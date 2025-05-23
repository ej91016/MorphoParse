# Changelog

All notable changes to this project will be documented in this file

## [v0.1.0b0] - 2025-05-07

### Note
- Start of formal version control

### Added
- Pre-release with support for **RAxML**, **RAxML-NG**, **IQ-TREE**, **PAUP\***, and **TNT**
- Robust parsing of morphological matrices:
  - Support for polymorphic encodings
  - Optional remapping of characters
  - Partition characters based on state space
- State-space-aware weighting using the `ln(r)` scheme