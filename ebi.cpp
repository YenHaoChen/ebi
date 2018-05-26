#include "ebi.h"

//#define NDEBUG // remove assert()

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <sstream>
#include <list>
using namespace std;

/*** base functions ******************************************************/

inline void ebi::base_initialization()
{ // initialize to zero
	sign = positive;
	N_bytes = 1;
	data = new uint8_t[1];
	data[0] = 0;
}

inline ebi ebi::base_addition(const ebi& a, const ebi& b) const
{ // add two positive number
	assert(a.sign == positive && b.sign == positive);

	if (a == 0)
		return b;
	if (b == 0)
		return a;

	bool result_sign = positive;
	unsigned int result_N_bytes = a.N_bytes>b.N_bytes? a.N_bytes : b.N_bytes;
	uint8_t result_data[MAX_DATA_LENGTH]; // we may get segmentation fault here

	int carry_in = 0;
	for (int i=0; i<(int)result_N_bytes; i++)
	{
		int sum = (int)a.get_byte(i) + (int)b.get_byte(i) + carry_in; // i may larger than the N_bytes of a and b
		result_data[i] = sum % 256;
		carry_in = sum / 256;
	}
	if (carry_in)
		result_data[result_N_bytes++] = carry_in;

	assert(result_N_bytes < MAX_DATA_LENGTH);
	return ebi(result_sign, result_N_bytes, result_data);
}

inline ebi ebi::base_subtraction(const ebi& a, const ebi& b) const
{ // subtracte big positive by small positive number
	assert(a.sign == positive && b.sign == positive);
	assert(a >= b);

	if (a == b)
		return 0;

	bool result_sign = positive;
	unsigned int result_N_bytes = a.N_bytes;
	uint8_t result_data[MAX_DATA_LENGTH]; // we may get segmentation fault here

	int borrow = 0;
	for (int i=0; i<(int)result_N_bytes; i++)
	{
		int sum = a.data[i] - b.get_byte(i) - borrow; // i may larger than b.N_bytes
		if (sum >= 0)
		{
			borrow = 0;
			result_data[i] = sum;
		}
		else
		{
			borrow = 1;
			result_data[i] = sum + 256;
		}
	}
	while (result_N_bytes > 1 && result_data[result_N_bytes-1] == 0)
		result_N_bytes--; // remove 0s, leave at least one digit

	return ebi(result_sign, result_N_bytes, result_data);
}

inline bool ebi::base_lessthan(const ebi& a, const ebi& b) const
{ // Compare two positive number
	assert(a.sign == positive && b.sign == positive);

	if (a.N_bytes < b.N_bytes)
		return true;
	if (a.N_bytes > b.N_bytes)
		return false;
	for (int i=a.N_bytes-1; i>=0; i--)
	{ // little endian
		if (a.data[i] < b.data[i])
			return true;
		if (a.data[i] > b.data[i])
			return false;
	}
	return false; // equal
}

/*** base functions *** end **********************************************/

//constructors

ebi::ebi()
{
	base_initialization(); // initialize to zero
}

ebi::ebi(int n)  //directly convert from an int
{
	sign = (n>=0 ? positive : negative);
	unsigned u = (unsigned)abs((long)n);
	N_bytes = (u&0xFF000000 ? 4 : (u&0x00FF0000 ? 3 : (u&0x0000FF00 ? 2 : 1)));
	data = new uint8_t[N_bytes];
	for (int i=0; i<(int)N_bytes; i++)
		data[i] = ((uint8_t*)(&u))[i];
}

ebi::ebi(bool isPositive, unsigned int nBlocks, uint8_t* rawData)
{
	if (nBlocks == 0 || (nBlocks==1 && rawData[0]==0))
	{
		base_initialization(); // initialize to zero
		return;
	}

	assert(nBlocks < MAX_DATA_LENGTH);
	sign = (isPositive ? positive : negative);
	N_bytes = nBlocks;
	while (N_bytes > 1 && rawData[N_bytes-1] == 0)
		N_bytes--; // remove 0s, leave at least one digit
	data = new uint8_t[N_bytes];
	for (int i=0; i<(int)N_bytes; i++)
		data[i] = rawData[i];
}

ebi::ebi(const ebi& bn)
{
	sign = bn.sign;
	N_bytes = bn.N_bytes;
	data = new uint8_t[N_bytes];
	for (int i=0; i<(int)N_bytes; i++)
		data[i] = bn.data[i];
}

ebi::ebi(const char* array)
{
	base_initialization(); // initialize to zero
	istringstream(array) >> *this;
}

//deconstructor

ebi::~ebi()
{
	assert(data);
	delete [] data;
	data = NULL;
}

//overloaded arithmetic operators as member functions
ebi ebi::operator+(const ebi &bn) const
{
	if (sign==negative && bn.sign==negative)
		return -base_addition(-(*this), -bn);

	if (sign==positive && bn.sign==negative && operator>=(-bn))
		return base_subtraction(*this, -bn);
	if (sign==positive && bn.sign==negative && operator<(-bn))
		return -base_subtraction(-bn, *this);

	if (sign==negative && bn.sign==positive && -(*this)<=bn)
		return base_subtraction(bn, -(*this));
	if (sign==negative && bn.sign==positive && -(*this)>bn)
		return -base_subtraction(-(*this), bn);

	return base_addition(*this, bn); // Adding two positive numbers
}

ebi ebi::operator-(const ebi &bn) const
{
	if (sign==negative && bn.sign==negative && -(*this)<=-bn)
		return base_subtraction(-bn, -(*this));
	if (sign==negative && bn.sign==negative && -(*this)>-bn)
		return -base_subtraction(-(*this), -bn);

	if (sign==positive && bn.sign==negative)
		return base_addition(*this, -bn);

	if (sign==negative && bn.sign==positive)
		return -base_addition(-(*this), bn);

	if (sign==positive && bn.sign==positive && operator<(bn))
		return -base_subtraction(bn, *this);
	//if (sign==positive && bn.sign==positive && operator>=(bn))
	return base_subtraction(*this, bn); // Subing positive numbers, big by small
}

ebi ebi::operator<<(unsigned int n) const
{
	if (n == 0 || operator==(0))
		return *this;

	bool result_sign = sign;
	unsigned int result_N_bytes = N_bytes + (n+BYTE_SIZE-1)/BYTE_SIZE;
	assert(result_N_bytes < MAX_DATA_LENGTH);
	uint8_t result_data[MAX_DATA_LENGTH]; // we may get segmentation fault here

	for (int i=(int)result_N_bytes-1; i>=(int)n/BYTE_SIZE; i--)
		result_data[i] = (get_byte(i-n/BYTE_SIZE)<<(n%BYTE_SIZE)) + (get_byte(i-(n+BYTE_SIZE-1)/BYTE_SIZE)>>(BYTE_SIZE-n%BYTE_SIZE));
	for (int i=0; i<(int)n/BYTE_SIZE; i++)
		result_data[i] = 0;

	return ebi(result_sign, result_N_bytes, result_data);
}

ebi ebi::operator>>(unsigned int n) const
{
	if (N_bytes <= n/BYTE_SIZE)
		return 0;

	bool result_sign = sign;
	unsigned int result_N_bytes = N_bytes - n/BYTE_SIZE;
	uint8_t result_data[MAX_DATA_LENGTH]; // we may get segmentation fault here
	for (int i=0; i<(int)result_N_bytes; i++)
		result_data[i] = (get_byte(i+n/BYTE_SIZE+1)<<(BYTE_SIZE-n%BYTE_SIZE)) + (get_byte(i+n/BYTE_SIZE)>>(n%BYTE_SIZE));

	return ebi(result_sign, result_N_bytes, (uint8_t*)result_data);
}

ebi ebi::operator*(const ebi &bn) const
{
	assert(N_bytes+bn.N_bytes < MAX_DATA_LENGTH);
	if (operator==(0) || bn==0)
		return 0;
	if (operator==(1))
		return bn;
	if (bn == 1)
		return *this;
	if (operator==(-1))
		return -bn;
	if (bn == -1)
		return -(*this);

	ebi result = 0;
	for (int i=0; i<(int)bn.N_bytes; i++)
	{
		ebi sum = 0;
		for (int j=0; j<(int)bn.data[i]; j++)
			sum += *this;
		result += sum << (i*BYTE_SIZE);
	}
	result.sign = (sign==bn.sign) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator/(const ebi &bn) const
{ //integer division: 3/2==1
	assert(bn != 0 && "divide by zero");
	if (bn == 1)
		return *this;

	ebi result = 0;
	ebi temp = abs(*this);
	ebi divisor = abs(bn);
	while (temp >= divisor)
	{
		unsigned int n = temp.N_bytes - divisor.N_bytes;
		ebi sub = divisor << (n*BYTE_SIZE);
		if (temp < sub)
		{
			assert(n != 0);
			n = n-1;
			sub = sub>>BYTE_SIZE;
		}

		temp -= sub;
		result += ebi(1)<<(n*BYTE_SIZE);
	}
	result.sign = (sign==bn.sign) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator%(const ebi &bn) const
{
	return (*this) - (((*this)/bn)*bn); // By C99, a == (a/b*b) + a%b
}

//interface functions

ostream& operator<<(ostream& out, const ebi& bn)
{
	if (bn.N_bytes==1 && bn.data[0]==0)
	{
		out << "0";
	}
	else
	{
		ios_base::fmtflags ff = out.flags();
		if (ff & out.dec)
		{
			if (bn.sign == ebi::negative)
				out << "-";
			list<int> dec_digits;
			for (ebi i=abs(bn); i!=0; i/=10)
				dec_digits.push_front((unsigned)(ebi(i%10).data[0]));
			for (list<int>::iterator it=dec_digits.begin(); it!=dec_digits.end(); it++)
				out << *it;
		}
		else if (ff & out.hex)
		{
			ios fmt_state(nullptr);
			fmt_state.copyfmt(out);
			ebi comp = bn; // output 2's complement
			if (comp.sign == ebi::negative)
			{ // translate to 2's complement
				comp.sign = ebi::positive; // currently not sure how to define this
				for (int i=0; i<(int)comp.N_bytes; i++)
					comp.data[i] = ~comp.data[i];
				comp = comp + 1; // should be considered together with sign
				if (!(comp.data[comp.N_bytes-1] & (1<<(BYTE_SIZE-1)))) // MSB is not 1, add another byte
					comp = comp + ((uint8_t(-1)) << (BYTE_SIZE*comp.N_bytes));
			}
			out << (unsigned)comp.data[comp.N_bytes-1]; // avoid 0 at beginning
			for (int i=comp.N_bytes-2; i>=0; i--) // little endian
				out << setw(BYTE_SIZE/4) << setfill('0') << (unsigned)comp.data[i];
			out.copyfmt(fmt_state);
		}
		else if (ff & out.oct)
		{ // this part is implemented, but never tested nor executed
			list<int> oct_digits;
			for (ebi i=abs(bn); i!=0; i=i>>3)
				oct_digits.push_front(i.data[0] % 8);
			for (list<int>::iterator it=oct_digits.begin(); it!=oct_digits.end(); it++)
				out << *it;
		}
		else
			assert(false && "Unknown fmtfl basefield");
	}
	return out;
}

istream& operator>>(istream& in, ebi& bn)
{
	char sign;
	while (isspace(sign=in.get()))
		continue; // read until first non-space character
	if (sign != '+' && sign != '-')
	{
		in.putback(sign);
		sign = '+';
	}

	char c_str[3] = {' ',' ','\0'};
	c_str[0] = in.get();
	c_str[1] = in.get();
	if (strcmp(c_str, "0x") && strcmp(c_str, "0X"))
	{ // decimal
		in.putback(c_str[1]);
		in.putback(c_str[0]);
		bn = 0;
		while(isdigit(c_str[0]=in.get()))
		{
			bn *= 10;
			bn += c_str[0] - '0';
		}
	}
	else
	{ // hexadecimal
		c_str[1] = '\0';
		list<uint8_t> hex_digits;
		while ((c_str[0]=in.get()) == '0') // read until first non-zero
			continue;
		in.putback(c_str[0]);
		while (isxdigit(c_str[0]=in.get())) // reverse required
			hex_digits.push_front((uint8_t)(strtol(c_str, NULL, 16)));
		assert(hex_digits.size() < MAX_DATA_LENGTH);

		if (hex_digits.size() == 0)
			bn = 0;
		else
		{
			bn.~ebi();
			if (hex_digits.size() % 2)
				hex_digits.push_back(0);
			bn.N_bytes = hex_digits.size() / 2;
			bn.data = new uint8_t[bn.N_bytes];
			int i = 0;
			list<uint8_t>::iterator it = hex_digits.begin();
			while (it != hex_digits.end())
			{
				bn.data[i] = *(it++);
				bn.data[i] = bn.data[i] + ((*(it++))<<4);
				i++;
			}
		}
	}

	bn.sign = ((bn==0 || sign=='+') ? ebi::positive : ebi::negative);

	return in;
}

bool ebi::operator==(const ebi &bn) const
{
	if (sign != bn.sign)
		return false;
	if (N_bytes != bn.N_bytes)
		return false;
	for (int i=N_bytes-1; i>=0; i--) // little endian
		if (data[i] != bn.data[i])
			return false;
	return true; // equal
}

bool ebi::operator!=(const ebi &bn) const
{
	return !operator==(bn);
}

bool ebi::operator<(const ebi &bn) const
{
	if (sign==negative && bn.sign==negative)
		return base_lessthan(-bn, -(*this));
	if (sign==positive && bn.sign==negative)
		return false;
	if (sign==negative && bn.sign==positive)
		return true;
	return base_lessthan(*this, bn); // Comparing two positive number
}

bool ebi::operator>(const ebi& bn) const
{
	return !operator<(bn) && !operator==(bn);
}

bool ebi::operator<=(const ebi &bn) const
{
	return operator<(bn) || operator==(bn);
}

bool ebi::operator>=(const ebi& bn) const
{
	return !operator<(bn);
}

unsigned int ebi::get_N_bytes() const
{
	return N_bytes;
}

uint8_t ebi::get_byte(unsigned int i) const
{ // also deal with i larger than N_bytes
	return i<N_bytes ? data[i] : 0;
}

ebi& ebi::operator=(const ebi& bn)
{
	this->~ebi();
	new (this) ebi(bn);
	return *this;
}

ebi& ebi::operator+=(const ebi& bn)
{
	return operator=( operator+(bn) );
}

ebi& ebi::operator-=(const ebi& bn)
{
	return operator=( operator-(bn) );
}

ebi& ebi::operator*=(const ebi& bn)
{
	return operator=( operator*(bn) );
}

ebi& ebi::operator/=(const ebi& bn)
{
	return operator=( operator/(bn) );
}

ebi ebi::operator-() const
{
	return ebi(!sign, N_bytes, data);
}

ebi& ebi::operator++()
{
	return operator+=(1);
}

ebi& ebi::operator--()
{
	return operator-=(1);
}

ebi ebi::operator++(int n)
{
	ebi temp = *this;
	operator+=(1);
	return temp;
}

ebi ebi::operator--(int n)
{
	ebi temp = *this;
	operator-=(1);
	return temp;
}

ebi abs(const ebi &bn)
{
	return bn>0? bn : -bn;
}

ebi pow(ebi base, unsigned exponent)
{ // base^exponent
	ebi n = 1;
	for (unsigned i=0; i<exponent; i++)
		n *= base;
	return n;
}

ebi operator+(int n, const ebi &bn)
{
	return ebi(n) + bn;
}

ebi operator-(int n, const ebi &bn)
{
	return ebi(n) - bn;
}

ebi operator*(int n, const ebi &bn)
{
	return ebi(n) * bn;
}

ebi operator/(int n, const ebi &bn)
{
	return ebi(n) / bn;
}

ebi operator%(int n, const ebi &bn)
{
	return ebi(n) % bn;
}

bool operator<(int n, const ebi &bn)
{
	return ebi(n) < bn;
}

bool operator!=(int n, const ebi &bn)
{
	return ebi(n) != bn;
}

bool operator==(int n, const ebi &bn)
{
	return ebi(n) == bn;
}

ebi::operator int() const
{
	assert(operator<(INT_MAX) && operator>(INT_MIN));
	int result = data[N_bytes-1];
	for (int i=N_bytes-2; i>=0; i--)
		result = (result<<BYTE_SIZE) + data[i];
	return sign==positive? result : -result;
}
