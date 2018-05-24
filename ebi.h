/**************************************************
  	Big-endian is implemented
	ebi(0) has positive sign and 1 digit
	int type is used for indexing, thus MAX_NUM_OF_BITS should not larger than INT_MAX
	I assume nobody will explicitly call the deconstructor
		i.e. should not exists a instance with NULL data array
**************************************************/

#ifndef __EXTREMELYEASYTOUSEBIGINTEGER__
#define __EXTREMELYEASYTOUSEBIGINTEGER__

#include <iostream>
using std::istream;
using std::ostream;
#include <cstdint>
#ifndef UINT8_MAX
typedef unsigned char uint8_t;
#endif

#define MAX_NUM_OF_BITS 10000

class ebi {
	public:
		//constructors and deconstructor
		ebi();
		ebi(int);  //directly convert from an int
		ebi(bool isPositive, unsigned int nXDigits, uint8_t* rawData);
		ebi(const ebi&);
		ebi(const char*);
		~ebi();

		//arithmetic operators
		ebi operator+(const ebi&) const;
		ebi operator-(const ebi&) const;
		ebi operator*(const ebi&) const;
		ebi operator/(const ebi&) const; //integer division: 3/2==1
		ebi operator%(const ebi&) const; //by C99, a == (a/b*b) + a%b
		ebi operator-() const;
		ebi operator<<(unsigned int) const;
		ebi operator>>(unsigned int) const;

		//interface functions
		bool get_sign() const; // discarded
		unsigned int get_N_xdigits() const;
		uint8_t get_data(unsigned int) const; // get_xdigit, operator[]
		friend ostream& operator<<(ostream&, const ebi&);
		friend istream& operator>>(istream&, ebi&);

		//comparison operators
		bool operator<(const ebi&) const;
		bool operator>(const ebi&) const;
		bool operator==(const ebi&) const;
		bool operator!=(const ebi&) const;
		bool operator<=(const ebi&) const;
		bool operator>=(const ebi&) const;

		//assignment operators
		ebi& operator=(const ebi&);
		ebi& operator+=(const ebi&);
		ebi& operator-=(const ebi&);
		ebi& operator*=(const ebi&);
		ebi& operator/=(const ebi&);

		//increment/decrement operators
		ebi& operator++();
		ebi& operator--();
		ebi operator++(int);
		ebi operator--(int);

		//cast
		explicit operator int() const;

	private:
		enum {negative, positive} sign;
		unsigned int N_xdigits;
		uint8_t *data;

		inline void base_initialization();
		inline ebi base_addition(const ebi &a, const ebi &b) const;
		inline ebi base_subtraction(const ebi& a, const ebi& b) const;
		inline bool base_lessthan(const ebi& a, const ebi& b) const;
};

ebi rand(const unsigned digits);

//math library
ebi abs(const ebi &);
ebi pow(ebi base, unsigned exponent);

ebi operator+(int n, const ebi &bn);
ebi operator-(int n, const ebi &bn);
ebi operator*(int n, const ebi &bn);
ebi operator/(int n, const ebi &bn);
ebi operator%(int n, const ebi &bn);
bool operator<(int n, const ebi &bn);
bool operator!=(int n, const ebi &bn);
bool operator==(int n, const ebi &bn);

#endif
