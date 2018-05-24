# An Extremely easy to use Big Integer (EBI) C++ library

EBI (Extremely easy to use Big Integer) is a multiple precision integer C++ library. The main feature is extremely easy to use. Source code is fairly simple. Basic C/C++ programming knowledge is enough to understand the implementation.

## Demo

C++ Code:
<pre><code>
	ebi a = "0xf1245ab3341ff3461818881767676819ee";
	ebi b;
	cout << "Please keyin a big number: ";
	cin >> b;
	cout << "Your number is " << b << endl;

	cout << "a+b = " << a+b << endl;
	cout << "a-b = " << a-b << endl;
	cout << "a*b = " << a*b << endl;
	cout << "a/b = " << a/b << endl;
	cout << "a\%b = " << a%b << endl;

	cout << "Counting from 0 to 9: ";
	for (ebi i=0; i<10; i++)
		cout << i << " ";
	cout << endl;
</code></pre>

Output:
<pre><code>
Please keyin a big number: 0x1234567
Your number is 0x1234567
a+b = 0xf1245ab3341ff3461818881767688b5f55
a-b = 0xf1245ab3341ff34618188817676644d487
a\*b = 0x1125db2ebc70b781e1a299501279ad4cf236994c2
a/b = 0xd3f0f41badb375e5b500492f1f04
a%b = 0x8a8b52
Counting from 0 to 9: 0x0 0x1 0x2 0x3 0x4 0x5 0x6 0x7 0x8 0x9
</code></pre>

## Features

- Support both decimal and hexadecimal format

## Implementation guidelines

- Other than constructors and internal functons (base\_\*), all of array pointers, uint8\_t \*data, are pointing allocated and independent memory location
