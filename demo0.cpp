#include <iostream>
using namespace std;

#include "ebi.h"

int main()
{
	ebi a = "0xf1245ab3341ff3461818881767676819ee";
	ebi b;
	cout << "Please keyin a big number: ";
	cin >> b;
	cout << "Your number is " << b << "(0x" << hex << b << ")" << endl;

	cout << "a+b = " << a+b << endl;
	cout << "a-b = " << a-b << endl;
	cout << "a*b = " << a*b << endl;
	cout << "a/b = " << a/b << endl;
	cout << "a\%b = " << a%b << endl;

	cout << "Counting from 0 to 9: ";
	for (ebi i=0; i<100; i++)
		cout << i << " ";
	cout << endl;

	return 0;
}
