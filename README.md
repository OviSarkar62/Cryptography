# Cryptography

An implementation of AES (Advanced Encryption Standard) built completely from scratch in Python — no cryptographic libraries (`cryptography`, `pycryptodome`, `hashlib`, etc.) used anywhere. Every lookup table, permutation, and round transformation is written out by hand so the internal mechanics of the cipher are fully visible.

## Contents

| File | Description |
|---|---|
| `AES_Encryption_Decryption.py` | AES-128 — 128-bit block, 128-bit key, 10 rounds. Includes the S-box/inverse S-box, key expansion (Rijndael schedule), SubBytes, ShiftRows, MixColumns, AddRoundKey, and their inverses. |

## How to Run

```bash
python3 AES_Encryption_Decryption.py
```

The script sets a sample plaintext and key at the top of the file, then prints:
- The plaintext and key (text/hex)
- The full key schedule (11 round keys)
- The state after every round of encryption
- The resulting ciphertext
- The state after every round of decryption
- The recovered plaintext, with a match check against the original

## Verification

- Validated against the official **NIST FIPS-197 Appendix C.1** test vector — the implementation reproduces the expected ciphertext byte-for-byte.
- Verified via round-trip: decrypting the generated ciphertext recovers the exact original plaintext.

## Purpose

Built as an educational project to understand block cipher internals — key scheduling, substitution-permutation networks, and Galois field arithmetic — by implementing every step manually rather than calling a library.

## Author

Ovi Sarkar
