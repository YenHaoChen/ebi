# An very Easy to use Big Integer (EBI) C++ library

The EBI (Easy to use Big Integer) is a multiple precision integer C++ library. The main feature is extremely easy to use. Source code is fairly simple. Basic C/C++ programming knowledge is enough to use and understand the implementation.

## Demo

C++ Code:
<pre><code>
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
	cout << "a%b = " << a%b << endl;
</code></pre>

Output:
<pre><code>
	Big number a = "12345678901234567890123456789"
	Please keyin a big number, e.g. 9876543210 or 0x24cb016ea:
	0x24cb016ea
	Your number, b, is 9876543210 (0x24cb016ea)

	a+b = 12345678901234567899999999999
	a-b = 12345678901234567880246913579
	a*b = 121932631124828532112482853211126352690
	a/b = 1249999988734374999
	a%b = 156249999
</code></pre>

## Features

- Source code is fairly simply, easy to use and understand
- Both decimal and hexadecimal formats are supported
- Small memory footprint

## Build

To try out the ebi C++ library, you only need to download three files, i.e. *ebi.h*, *ebi.cpp* and *demo0.cpp*.

Under a Linux environment:
<pre><code>
	$ gcc -std=c++11 -o demo demo0.cpp ebi.cpp
	$ ./demo
</code></pre>

Under a Windows environment with CodeBlocks 17.12:
<pre><code>
	Settings -> Compiler...
		Select "Compiler settings"
		Select "Compiler Flags"
		Select "Have g++ follow the C++11 ISO C++ language standard [-std=c++11]"
		Press OK buttom
	Build -> Build and run (F9)
</code></pre>

---------------------------------------------------------------------

## Implementation details

- The c++11 features are used in the ebi library which has been supported by most of widely used compilers ([C++ compiler support](https://en.cppreference.com/w/cpp/compiler_support "C++ compiler support"))
- Little-endian is implemented, i.e. highest address as the Most Significant Byte (MSB)
- The backbone data structure is a uint8\_t array (or unsigned char array); thus, the ebi library has following advantages:
	- The ebi library is efficient in terms of memory usage
	- Shift operators have the same physical memory meaning as int type

## Implementation guidelines

- Other than constructors and internal functons (base\_\*), all of array pointers, uint8\_t \*data, are pointing allocated and independent memory location

## Development platform

- CentOS Linux 7
- x86-64 ISA
- gcc (GCC) 4.8.5

## Issue raising is welcomed

The first motivation of ebi C++ library is simple, easy to use and understand. And, yea, issue raising is welcomed. ^^
