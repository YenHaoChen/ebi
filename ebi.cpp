// By vegetablebird 2018.03.09

#include "ebi.h"

//#define NDEBUG // remove assert()
//#define DEBUG_BIGNUMBER

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <sstream>
using namespace std;

//constructors

#define INITIALIZE() { \
	sgn = positive; \
	num_of_bits = 1; \
	data = new uint8_t[1]; \
	data[0] = 0; \
}


ebi::ebi()
{
	INITIALIZE(); // initialize to zero
}

ebi::ebi(int n)  //directly convert from an int
{
	assert(abs(n) >= 0 && "Notice: abs(INT_MIN) is negative, currently only support to INT_MIN+1");
#ifdef DEBUG_BIGNUMBER
	cout << "Creating (int)" << n << endl;
#endif
	if (n == 0)
	{
		INITIALIZE(); // initialize to zero
		return;
	}

	sgn = n>=0? positive : negative;
	n = abs(n);
	double num_of_bits_fp = ceil(log2(n)/4.0);
	assert(num_of_bits_fp == num_of_bits_fp);
	num_of_bits = (unsigned int)num_of_bits_fp;
	if (1<<(num_of_bits*4) == n) // power of 16
	{ // n=1, 16, 16^2, 16^3, ...
		num_of_bits++;
		assert(num_of_bits < MAX_NUM_OF_BITS);
		data = new uint8_t[num_of_bits];
		for (int i=0; i<(int)num_of_bits-1; i++)
			data[i] = 0;
		data[num_of_bits-1] = 1;
	}
	else
	{
		assert(num_of_bits < MAX_NUM_OF_BITS);
		data = new uint8_t[num_of_bits];
		for (int i=0; i<(int)num_of_bits; i++)
		{ // big endian
			data[i] = n % 16;
			n = n / 16;
		}
	}
}

ebi::ebi(bool s, unsigned int n, uint8_t* d)
{
	if (n == 0 || (n==1 && d[0]==0))
	{
		INITIALIZE(); // initialize to zero
		return;
	}

	assert(n < MAX_NUM_OF_BITS);
	sgn = s;
	num_of_bits = n;
	data = new uint8_t[num_of_bits];
	while (num_of_bits > 1 && d[num_of_bits-1] == 0)
		num_of_bits--; // remove 0s, leave at least one digit
	for (int i=0; i<(int)num_of_bits; i++)
		data[i] = d[i];
}

ebi::ebi(const ebi& bn)
{
#ifdef DEBUG_BIGNUMBER
	cout << "Copying " << bn << endl;
#endif
	sgn = bn.sgn;
	num_of_bits = bn.num_of_bits;
	data = new uint8_t[num_of_bits];
	for (int i=0; i<(int)num_of_bits; i++)
		data[i] = bn.data[i];
}

ebi::ebi(const char* array)
{
	INITIALIZE(); // initialize to zero
	istringstream(array) >> *this;
}

//deconstructor

ebi::~ebi()
{
#ifdef DEBUG_BIGNUMBER
	cout << "Deleting " << *this << endl;
#endif
	assert(data);
	delete [] data;
	data = NULL;
}

//overloaded arithmetic operators as member functions
ebi ebi::operator+(ebi bn)
{ // I hope compiler can optimize this...
#ifdef DEBUG_BIGNUMBER
	cout << "Adding " << *this << " " << bn << endl;
#endif
	if (isNegative() && bn.isNegative())
		return -base_addition(-(*this), -bn);

	if (isPositive() && bn.isNegative() && operator>=(-bn))
		return base_subtraction(*this, -bn);
	if (isPositive() && bn.isNegative() && operator<(-bn))
		return -base_subtraction(-bn, *this);

	if (isNegative() && bn.isPositive() && -(*this)<=bn)
		return base_subtraction(bn, -(*this));
	if (isNegative() && bn.isPositive() && -(*this)>bn)
		return -base_subtraction(-(*this), bn);

	return base_addition(*this, bn); // Adding two positive numbers
}

ebi ebi::operator-(ebi bn)
{ // I hope compiler can optimize this...
#ifdef DEBUG_BIGNUMBER
	cout << "Subtracting " << *this << " " << bn << endl;
#endif
	if (isNegative() && bn.isNegative() && -(*this)<=-bn)
		return base_subtraction(-bn, -(*this));
	if (isNegative() && bn.isNegative() && -(*this)>-bn)
		return -base_subtraction(-(*this), -bn);

	if (isPositive() && bn.isNegative())
		return base_addition(*this, -bn);

	if (isNegative() && bn.isPositive())
		return -base_addition(-(*this), bn);

	if (isPositive() && bn.isPositive() && operator<(bn))
		return -base_subtraction(bn, *this);
	//if (isPositive() && bn.isPositive() && operator>=(bn))
	return base_subtraction(*this, bn); // Subing positive numbers, big by small
}

ebi ebi::operator<<(unsigned int n) const
{
	assert (n % 4 == 0); // FIXME
	n /= 4;
#ifdef DEBUG_BIGNUMBER
	cout << "Shifting left " << *this << " by " << n << endl;
#endif
	assert(n < MAX_NUM_OF_BITS);
	if (n == 0 || operator==(0))
		return *this;

	bool result_sgn = sgn;
	unsigned int result_num_of_bits = num_of_bits + n;
	assert(result_num_of_bits < MAX_NUM_OF_BITS);
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	for (int i=num_of_bits-1; i>=0; i--)
		result_data[i+n] = data[i];
	for (int i=0; i<(int)n; i++)
		result_data[i] = 0;

	return ebi(result_sgn, result_num_of_bits, result_data);
}

ebi ebi::operator>>(unsigned int n) const
{
	assert (n % 4 == 0); // FIXME
	n /= 4;
#ifdef DEBUG_BIGNUMBER
	cout << "Shifting right " << *this << " by " << n << endl;
#endif
	if (num_of_bits <= n)
		return 0;

	bool result_sgn = sgn;
	unsigned int result_num_of_bits = num_of_bits - n;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	for (int i=n; i<(int)num_of_bits; i++)
		result_data[i-n] = data[i];

	return ebi(result_sgn, result_num_of_bits, result_data);
}

ebi ebi::operator*(ebi bn)
{
#ifdef DEBUG_BIGNUMBER
	cout << "Multiplying " << *this << " " << bn << endl;
#endif
	assert(num_of_bits+bn.num_of_bits < MAX_NUM_OF_BITS);
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
	for (int i=0; i<(int)bn.num_of_bits; i++)
	{
		ebi sum = 0;
		for (int j=0; j<(int)bn.data[i]; j++)
			sum += *this;
		result += sum << (i*4);
	}
	result.sgn = (sgn==bn.sgn) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator/(ebi bn) //integer division: 3/2==1
{
#ifdef DEBUG_BIGNUMBER
	cout << "Dividing " << *this << " " << bn << endl;
#endif
	assert(bn != 0 && "divide by zero");
	if (bn == 1)
		return *this;

	ebi result = 0;
	ebi temp = abs(*this);
	ebi divisor = abs(bn);
	while (temp >= divisor)
	{
		unsigned int n = temp.num_of_bits - divisor.num_of_bits;
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
	result.sgn = (sgn==bn.sgn) ? positive : negative; // must be assigned lastly
	return result;
}

ebi ebi::operator%(ebi bn)
{ // not support modulus by a negative
#ifdef DEBUG_BIGNUMBER
	cout << "Modular " << *this << " " << bn << endl;
#endif
	assert(bn > 0 && "Modulo by non-positive number");
	if (bn == 1)
		return 0;

	ebi result = abs(*this);
	while (result >= bn)
	{
		unsigned int n = result.num_of_bits - bn.num_of_bits;
		ebi sub = bn << (n*4);
		result -= (result>=sub) ? sub : sub>>4;
	}
	if (isNegative())
		result = bn - result;
	return result;
}

//interface functions

void ebi::Print()
{
	cout << *this;
}

void ebi::GetData(bool& s, unsigned int& n, uint8_t* d)
{
	s = sgn;
	n = num_of_bits;
	d = data;
}

ostream& operator<<(ostream& out, const ebi& bn)
{
	if (bn.sgn == ebi::negative)
		out << "-";
	out << "0x";
	for (int i=bn.num_of_bits-1; i>=0; i--) // big endian
		out << hex << (int)bn.data[i];
	return out << dec;
}

istream& operator>>(istream& in, ebi& bn)
{
	bool sgn;
	char c;
	do {
		c = in.get();
//	} while (c==' ' || c=='\n'); // read until first non-space character
	} while (c!='+' && c!='-' && !isxdigit(c)); // read until first valid character
	if (c == '+')
		sgn = ebi::positive;
	else if (c == '-')
		sgn = ebi::negative;
	else if (isxdigit(c))
	{
		sgn = ebi::positive;
		in.putback(c);
	}
	else
	{ // invalid character, intinal to 0
		cout << "invalid character: '" << c << "', inital to 0" << endl;
		bn = 0;
		return in;
	}

	unsigned int num_of_bits = 0;
	uint8_t data[MAX_NUM_OF_BITS]; // we may get segmentation fault here
	char c_str[2] = {' ','\0'};
	c_str[0] = in.get();
	c_str[1] = in.get();
	assert((!strncmp(c_str, "0x", 2) || !strncmp(c_str, "0X", 2)) && "Only support hex string input, e.g. 0x333");
//	in.putback(c_str[1]);
//	in.putback(c_str[2]);
	while (isxdigit(c_str[0]=in.get()))
	{ // input is little endian, need to turn into big endian
		data[num_of_bits++] = (uint8_t)strtol(c_str, NULL, 16);
		assert(num_of_bits < MAX_NUM_OF_BITS);
	}
	for (int i=0; i<(int)num_of_bits/2; i++)
	{ // input is little endian, swap into big endian
		uint8_t temp = data[i];
		data[i] = data[num_of_bits-i-1];
		data[num_of_bits-i-1] = temp;
	}

	bn = ebi(sgn, num_of_bits, data);
	return in;
}

inline bool ebi::isPositive() const
{
	return sgn==positive;
}

inline bool ebi::isNegative() const
{
	return sgn==negative;
}

bool ebi::operator==(const ebi &bn) const
{
	if (sgn != bn.sgn)
		return false;
	if (num_of_bits != bn.num_of_bits)
		return false;
	for (int i=num_of_bits-1; i>=0; i--) // big endian
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
	if (sgn==negative && bn.sgn==negative)
		return base_lessthan(-bn, -(*this));
	if (sgn==positive && bn.sgn==negative)
		return false;
	if (sgn==negative && bn.sgn==positive)
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

bool ebi::get_sgn() const
{
	return sgn;
}

unsigned int ebi::get_num_of_bits() const
{
	return num_of_bits;
}

uint8_t ebi::get_data(unsigned int i) const
{ // also deal with i larger than num_of_bits
	return i<num_of_bits ? data[i] : 0;
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
#ifdef DEBUG_BIGNUMBER
	cout << "Negating " << (*this) << " -> " << ebi(!sgn, num_of_bits, data) << endl;
#endif
	return ebi(!sgn, num_of_bits, data);
}

inline ebi ebi::base_addition(const ebi& a, const ebi& b) const
{ // add two positive number
	assert(a.isPositive());
	assert(b.isPositive());

	if (a == 0)
		return b;
	if (b == 0)
		return a;

	bool result_sgn = positive;
	unsigned int result_num_of_bits = a.num_of_bits>b.num_of_bits? a.num_of_bits : b.num_of_bits;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	int carry_in = 0;
	for (int i=0; i<(int)result_num_of_bits; i++)
	{
		int sum = a.get_data(i) + b.get_data(i) + carry_in; // i may larger than the num_of_bits of a and b
		result_data[i] = sum % 16;
		carry_in = sum / 16;
	}
	if (carry_in)
		result_data[result_num_of_bits++] = carry_in;

	assert(result_num_of_bits < MAX_NUM_OF_BITS);
	return ebi(result_sgn, result_num_of_bits, result_data);
}

inline ebi ebi::base_subtraction(const ebi& a, const ebi& b) const
{ // subtracte big positive by small positive number
	assert(a.isPositive());
	assert(b.isPositive());
	assert(a >= b);

	if (a == b)
		return 0;

	bool result_sgn = positive;
	unsigned int result_num_of_bits = a.num_of_bits;
	uint8_t result_data[MAX_NUM_OF_BITS]; // we may get segmentation fault here

	int borrow = 0;
	for (int i=0; i<(int)result_num_of_bits; i++)
	{
		int sum = a.data[i] - b.get_data(i) - borrow; // i may larger than b.num_of_bits
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
	while (result_num_of_bits > 1 && result_data[result_num_of_bits-1] == 0)
		result_num_of_bits--; // remove 0s, leave at least one digit

	return ebi(result_sgn, result_num_of_bits, result_data);
}

inline bool ebi::base_lessthan(const ebi& a, const ebi& b) const
{ // Compare two positive number
	assert(a.isPositive());
	assert(b.isPositive());
#ifdef DEBUG_BIGNUMBER
	cout << "Comparing " << a << " < " << b << " ? " << endl;
#endif

	if (a.num_of_bits < b.num_of_bits)
		return true;
	if (a.num_of_bits > b.num_of_bits)
		return false;
	for (int i=a.num_of_bits-1; i>=0; i--)
	{ // big endian
		if (a.data[i] < b.data[i])
			return true;
		if (a.data[i] > b.data[i])
			return false;
	}
	return false; // equal
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

ebi operator+(int n, ebi bn)
{
	return ebi(n) + bn;
}

ebi operator-(int n, ebi bn)
{
	return ebi(n) - bn;
}

ebi operator*(int n, ebi bn)
{
	return ebi(n) * bn;
}

ebi operator/(int n, ebi bn)
{
	return ebi(n) / bn;
}

ebi operator%(int n, ebi bn)
{
	return ebi(n) % bn;
}

bool operator<(int n, ebi bn)
{
	return ebi(n) < bn;
}

bool operator!=(int n, ebi bn)
{
	return ebi(n) != bn;
}

bool operator==(int n, ebi bn)
{
	return ebi(n) == bn;
}

ebi::operator int() const
{
	assert(operator<(INT_MAX) && operator>(INT_MIN+1));
	int result = data[num_of_bits-1];
	for (int i=num_of_bits-2; i>=0; i--)
		result = (result<<4) + data[i];
	return sgn==positive? result : -result;
}
