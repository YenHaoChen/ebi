// By vegetablebird 2018.03.09

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
	N_xdigits = 1;
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
	unsigned int result_N_xdigits = a.N_xdigits>b.N_xdigits? a.N_xdigits : b.N_xdigits;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	int carry_in = 0;
	for (int i=0; i<(int)result_N_xdigits; i++)
	{
		int sum = a.get_data(i) + b.get_data(i) + carry_in; // i may larger than the N_xdigits of a and b
		result_data[i] = sum % 16;
		carry_in = sum / 16;
	}
	if (carry_in)
		result_data[result_N_xdigits++] = carry_in;

	assert(result_N_xdigits < MAX_NUM_OF_BITS);
	return ebi(result_sign, result_N_xdigits, result_data);
}

inline ebi ebi::base_subtraction(const ebi& a, const ebi& b) const
{ // subtracte big positive by small positive number
	assert(a.sign == positive && b.sign == positive);
	assert(a >= b);

	if (a == b)
		return 0;

	bool result_sign = positive;
	unsigned int result_N_xdigits = a.N_xdigits;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	int borrow = 0;
	for (int i=0; i<(int)result_N_xdigits; i++)
	{
		int sum = a.data[i] - b.get_data(i) - borrow; // i may larger than b.N_xdigits
		if (sum >= 0)
		{
			borrow = 0;
			result_data[i] = sum;
		}
		else
		{
			borrow = 1;
			result_data[i] = sum + 16;
		}
	}
	while (result_N_xdigits > 1 && result_data[result_N_xdigits-1] == 0)
		result_N_xdigits--; // remove 0s, leave at least one digit

	return ebi(result_sign, result_N_xdigits, result_data);
}

inline bool ebi::base_lessthan(const ebi& a, const ebi& b) const
{ // Compare two positive number
	assert(a.sign == positive && b.sign == positive);

	if (a.N_xdigits < b.N_xdigits)
		return true;
	if (a.N_xdigits > b.N_xdigits)
		return false;
	for (int i=a.N_xdigits-1; i>=0; i--)
	{ // big endian
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
	assert(abs(n) >= 0 && "Notice: abs(INT_MIN) is negative, currently only support to INT_MIN+1");
	if (n == 0)
	{
		base_initialization(); // initialize to zero
		return;
	}

	sign = n>=0? positive : negative;
	n = abs(n);
	double N_xdigits_fp = ceil(log2(n)/4.0);
	assert(N_xdigits_fp == N_xdigits_fp);
	N_xdigits = (unsigned int)N_xdigits_fp;
	if (1<<(N_xdigits*4) == n) // power of 16
	{ // n=1, 16, 16^2, 16^3, ...
		N_xdigits++;
		assert(N_xdigits < MAX_NUM_OF_BITS);
		data = new uint8_t[N_xdigits];
		for (int i=0; i<(int)N_xdigits-1; i++)
			data[i] = 0;
		data[N_xdigits-1] = 1;
	}
	else
	{
		assert(N_xdigits < MAX_NUM_OF_BITS);
		data = new uint8_t[N_xdigits];
		for (int i=0; i<(int)N_xdigits; i++)
		{ // big endian
			data[i] = n % 16;
			n = n / 16;
		}
	}
}

ebi::ebi(bool isPositive, unsigned int nXDigits, uint8_t* rawData)
{
	if (nXDigits == 0 || (nXDigits==1 && rawData[0]==0))
	{
		base_initialization(); // initialize to zero
		return;
	}

	assert(nXDigits < MAX_NUM_OF_BITS);
	sign = (isPositive ? positive : negative);
	N_xdigits = nXDigits;
	data = new uint8_t[N_xdigits];
	while (N_xdigits > 1 && rawData[N_xdigits-1] == 0)
		N_xdigits--; // remove 0s, leave at least one digit
	for (int i=0; i<(int)N_xdigits; i++)
		data[i] = rawData[i];
}

ebi::ebi(const ebi& bn)
{
	sign = bn.sign;
	N_xdigits = bn.N_xdigits;
	data = new uint8_t[N_xdigits];
	for (int i=0; i<(int)N_xdigits; i++)
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
{ // I hope compiler can optimize this...
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
{ // I hope compiler can optimize this...
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
	assert (n % 4 == 0); // FIXME
	n /= 4;
	assert(n < MAX_NUM_OF_BITS);
	if (n == 0 || operator==(0))
		return *this;

	bool result_sign = sign;
	unsigned int result_N_xdigits = N_xdigits + n;
	assert(result_N_xdigits < MAX_NUM_OF_BITS);
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	for (int i=N_xdigits-1; i>=0; i--)
		result_data[i+n] = data[i];
	for (int i=0; i<(int)n; i++)
		result_data[i] = 0;

	return ebi(result_sign, result_N_xdigits, result_data);
}

ebi ebi::operator>>(unsigned int n) const
{
	assert (n % 4 == 0); // FIXME
	n /= 4;
	if (N_xdigits <= n)
		return 0;

	bool result_sign = sign;
	unsigned int result_N_xdigits = N_xdigits - n;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	for (int i=n; i<(int)N_xdigits; i++)
		result_data[i-n] = data[i];

	return ebi(result_sign, result_N_xdigits, result_data);
}

ebi ebi::operator*(const ebi &bn) const
{
	assert(N_xdigits+bn.N_xdigits < MAX_NUM_OF_BITS);
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
	for (int i=0; i<(int)bn.N_xdigits; i++)
	{
		ebi sum = 0;
		for (int j=0; j<(int)bn.data[i]; j++)
			sum += *this;
		result += sum << (i*4);
	}
	result.sign = (sign==bn.sign) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator/(const ebi &bn) const //integer division: 3/2==1
{
	assert(bn != 0 && "divide by zero");
	if (bn == 1)
		return *this;

	ebi result = 0;
	ebi temp = abs(*this);
	ebi divisor = abs(bn);
	while (temp >= divisor)
	{
		unsigned int n = temp.N_xdigits - divisor.N_xdigits;
		ebi sub = divisor << (n*4);
		if (temp >= sub)
		{
			temp -= sub;
			result += ebi(1)<<(n*4);
		}
		else
		{
			temp -= sub>>4;
			result += ebi(1)<<((n-1)*4);
		}
	}
	result.sign = (sign==bn.sign) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator%(const ebi &bn) const
{
	return (*this) - (((*this)/bn)*bn); // By C99, a == (a/b*b) + a%b
}

//interface functions

void ebi::Print()
{
	cout << *this;
}

void ebi::GetData(bool& s, unsigned int& n, uint8_t* d)
{
	s = sign;
	n = N_xdigits;
	d = data;
}

ostream& operator<<(ostream& out, const ebi& bn)
{
	ios_base::fmtflags ff = out.flags();
	if (bn.sign == ebi::negative)
		out << "-";
	if (ff & out.dec)
	{
		list<int> dec_digits;
		for (ebi i=abs(bn); i!=0; i/=10)
			dec_digits.push_front((unsigned)(ebi(i%10).data[0]));
		for (list<int>::iterator it=dec_digits.begin(); it!=dec_digits.end(); it++)
			cout << *it;
	}
	else if (ff & out.hex)
	{
		for (int i=bn.N_xdigits-1; i>=0; i--) // big endian
			out << hex << (unsigned)bn.data[i];
	}
	else if (ff & out.oct)
		assert(false && "Error: octadecimal output format not supported");
	else
		assert(false && "Unknown fmtfl basefield");
	return out << dec;
}

istream& operator>>(istream& in, ebi& bn)
{
	bool sign;
	char c;
	do {
		c = in.get();
//	} while (c==' ' || c=='\n'); // read until first non-space character
	} while (c!='+' && c!='-' && !isxdigit(c)); // read until first valid character
	if (c == '+')
		sign = ebi::positive;
	else if (c == '-')
		sign = ebi::negative;
	else if (isxdigit(c))
	{
		sign = ebi::positive;
		in.putback(c);
	}
	else
	{ // invalid character, intinal to 0
		cout << "invalid character: '" << c << "', inital to 0" << endl;
		bn = 0;
		return in;
	}

	unsigned int N_xdigits = 0;
	uint8_t data[MAX_NUM_OF_BITS]; // we may get segmentation fault here
	char c_str[2] = {' ','\0'};
	c_str[0] = in.get();
	c_str[1] = in.get();
	assert((!strncmp(c_str, "0x", 2) || !strncmp(c_str, "0X", 2)) && "Only support hex string input, e.g. 0x333"); // FIXME
//	in.putback(c_str[1]);
//	in.putback(c_str[2]);
	while (isxdigit(c_str[0]=in.get()))
	{ // input is little endian, need to turn into big endian
		data[N_xdigits++] = (uint8_t)strtol(c_str, NULL, 16);
		assert(N_xdigits < MAX_NUM_OF_BITS);
	}
	for (int i=0; i<(int)N_xdigits/2; i++)
	{ // input is little endian, swap into big endian
		uint8_t temp = data[i];
		data[i] = data[N_xdigits-i-1];
		data[N_xdigits-i-1] = temp;
	}

	bn = ebi(sign, N_xdigits, data);
	return in;
}

bool ebi::operator==(const ebi &bn) const
{
	if (sign != bn.sign)
		return false;
	if (N_xdigits != bn.N_xdigits)
		return false;
	for (int i=N_xdigits-1; i>=0; i--) // big endian
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

bool ebi::get_sign() const
{
	return sign;
}

unsigned int ebi::get_N_xdigits() const
{
	return N_xdigits;
}

uint8_t ebi::get_data(unsigned int i) const
{ // also deal with i larger than N_xdigits
	return i<N_xdigits ? data[i] : 0;
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
	return ebi(!sign, N_xdigits, data);
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

ebi rand(const unsigned digits)
{
	ebi n = 0;
	for (unsigned i=0; i<digits; i++)
		n = (n<<4) + (rand()%16);
	return n;
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
	assert(operator<(INT_MAX) && operator>(INT_MIN+1));
	int result = data[N_xdigits-1];
	for (int i=N_xdigits-2; i>=0; i--)
		result = (result<<4) + data[i];
	return sign==positive? result : -result;
}
