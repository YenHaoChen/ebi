# Changelog

## [Unreleased]

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
