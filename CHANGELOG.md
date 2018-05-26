# Changelog

## [Unreleased]

## 2018.05.26
- 2''s complement hexadecimal output format
- Bug fixing
- Regression framework is completed; thus, this should be the first stable version

## 2018.05.25
- sed -i 's/MAX_NUM_OF_BITS/MAX_DATA_LENGTH/g' *
- sed -i 's/N_xdigits/N_bytes/g' *
- Reduce memory footprint by half
- Deal with corner case, INT_MIN
- Rewrite ebi(int)
- sed -i 's/get_byte/get_byte/g' *
- Add MACRO: #define EasyBigInteger ebi
- Update README.md
- sed -i 's/block/byte/g' *
- Complete shift operators, operator<<(int) and operator>>(int)
- Implement octadeciaml output, but not tested nor executed
- Add macro BYTE_SIZE

## 2018.05.24
- Check uint8_t exists by #ifndef UINT8_MAX
- Specify using std::istream in ebi.h
- Specify using std::ostream in ebi.h
- Add implementation guideline to README.md
- Support deciaml and hexadecimal input/output format

## 2018.05.24
- Add Demo and Features in README.md
- Modify Operator%() by C99 definition, a == (a/b*b) + a%b
- Rename num_of_bits to N_xdigits
- Rename __BIGNUMBER__ to __EXTREMELYEASYTOUSEBIGINTEGER__
- Remove DEBUG_BIGNUMBER
- Support decimal and hexadecimal output
- Remove isPositive() and isNegative()
- Rename bool sgn to enum sign
- Remove constants positive and negative
- Modify INITIALIZE() to base_initialization()
