/**************************************************
  	Little-endian is implemented
	ebi(0) has positive sign and 1 digit
	int type is used for indexing, thus MAX_DATA_LENGTH should not larger than INT_MAX
	I assume nobody will explicitly call the deconstructor
		i.e. each instance has a valid data array
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

#define EasyBigInteger ebi

#define MAX_DATA_LENGTH 10000
#define BYTE_SIZE 8

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
		unsigned int get_N_bytes() const;
		uint8_t get_byte(unsigned int) const;
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
		unsigned int N_bytes;
		uint8_t *data;

		inline void base_initialization();
		inline ebi base_addition(const ebi &a, const ebi &b) const;
		inline ebi base_subtraction(const ebi& a, const ebi& b) const;
		inline bool base_lessthan(const ebi& a, const ebi& b) const;
};

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
