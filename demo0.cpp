#include <iostream>
using namespace std;

#include "ebi.h"

int main()
{
	ebi a = "12345678901234567890123456789";
	cout << "Big number a = \"" << a << "\"" << endl;
	ebi b;
	cout << "Please keyin a big number, e.g. 9876543210 or 0x24cb016ea:" << endl;
	cin >> b;
	cout << "Your number, b, is " << b << " (0x" << hex << b << ")" << endl;
	cout << endl;
	cout << "a+b = " << a+b << endl;
	cout << "a-b = " << a-b << endl;
	cout << "a*b = " << a*b << endl;
	cout << "a/b = " << a/b << endl;
	cout << "a\%b = " << a%b << endl;

	return 0;
}
