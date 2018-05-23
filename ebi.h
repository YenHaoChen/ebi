//Big number class definition. This is just an example.

/**************************************************
By vegetablebird 2018.03.09
  	Big-endian is implemented
	ebi(0) has positive sign and 1 digit
	int type is used for indexing, thus MAX_NUM_OF_BITS should not larger than INT_MAX
	I assume nobody will explicitly call the deconstructure
		i.e. should not exists a instance with NULL data array
**************************************************/

#ifndef __BIGNUMBER__
#define __BIGNUMBER__

#include <iostream>
//#include <cstdint>
using namespace std;

typedef unsigned char uint8_t;

#define MAX_NUM_OF_BITS 10000

class ebi {
	private:
		bool sgn; // sign
		unsigned int num_of_bits; // N_xdigits
		uint8_t *data;

		const static bool positive = true;
		const static bool negative = false;

		inline ebi base_addition(const ebi &a, const ebi &b) const;
		inline ebi base_subtraction(const ebi& a, const ebi& b) const;
		inline bool base_lessthan(const ebi& a, const ebi& b) const;

		inline bool isPositive() const;
		inline bool isNegative() const;

	public:
		//constructors
		ebi();
		ebi(int);  //directly convert from an int
		ebi(bool, unsigned int, uint8_t*);
		ebi(const ebi&);
		ebi(const char*);

		//deconstructor
		~ebi();

		//overloaded arithmetic operators as member functions
		ebi operator+(ebi);
		ebi operator-(ebi);
		ebi operator*(ebi);
		ebi operator/(ebi); //integer division: 3/2==1
		ebi operator%(ebi);
		ebi operator-() const;
		ebi operator<<(unsigned int) const;
		ebi operator>>(unsigned int) const;

		//interface functions
		void Print();
		void GetData(bool&, unsigned int&, uint8_t*);
		bool get_sgn() const;
		unsigned int get_num_of_bits() const;
		uint8_t get_data(unsigned int) const;
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
};

ebi rand(const unsigned digits);

//math library
ebi abs(const ebi &);
ebi pow(ebi base, unsigned exponent);

ebi operator+(int n, ebi bn);
ebi operator-(int n, ebi bn);
ebi operator*(int n, ebi bn);
ebi operator/(int n, ebi bn);
ebi operator%(int n, ebi bn);
bool operator<(int n, ebi bn);
bool operator!=(int n, ebi bn);
bool operator==(int n, ebi bn);

#endif
