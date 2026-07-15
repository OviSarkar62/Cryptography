# ---------------------------------------------------------------------------
# AES S-box and its inverse (the fixed substitution tables from FIPS-197)
# ---------------------------------------------------------------------------
S_BOX = [
    0x63,0x7c,0x77,0x7b,0xf2,0x6b,0x6f,0xc5,0x30,0x01,0x67,0x2b,0xfe,0xd7,0xab,0x76,
    0xca,0x82,0xc9,0x7d,0xfa,0x59,0x47,0xf0,0xad,0xd4,0xa2,0xaf,0x9c,0xa4,0x72,0xc0,
    0xb7,0xfd,0x93,0x26,0x36,0x3f,0xf7,0xcc,0x34,0xa5,0xe5,0xf1,0x71,0xd8,0x31,0x15,
    0x04,0xc7,0x23,0xc3,0x18,0x96,0x05,0x9a,0x07,0x12,0x80,0xe2,0xeb,0x27,0xb2,0x75,
    0x09,0x83,0x2c,0x1a,0x1b,0x6e,0x5a,0xa0,0x52,0x3b,0xd6,0xb3,0x29,0xe3,0x2f,0x84,
    0x53,0xd1,0x00,0xed,0x20,0xfc,0xb1,0x5b,0x6a,0xcb,0xbe,0x39,0x4a,0x4c,0x58,0xcf,
    0xd0,0xef,0xaa,0xfb,0x43,0x4d,0x33,0x85,0x45,0xf9,0x02,0x7f,0x50,0x3c,0x9f,0xa8,
    0x51,0xa3,0x40,0x8f,0x92,0x9d,0x38,0xf5,0xbc,0xb6,0xda,0x21,0x10,0xff,0xf3,0xd2,
    0xcd,0x0c,0x13,0xec,0x5f,0x97,0x44,0x17,0xc4,0xa7,0x7e,0x3d,0x64,0x5d,0x19,0x73,
    0x60,0x81,0x4f,0xdc,0x22,0x2a,0x90,0x88,0x46,0xee,0xb8,0x14,0xde,0x5e,0x0b,0xdb,
    0xe0,0x32,0x3a,0x0a,0x49,0x06,0x24,0x5c,0xc2,0xd3,0xac,0x62,0x91,0x95,0xe4,0x79,
    0xe7,0xc8,0x37,0x6d,0x8d,0xd5,0x4e,0xa9,0x6c,0x56,0xf4,0xea,0x65,0x7a,0xae,0x08,
    0xba,0x78,0x25,0x2e,0x1c,0xa6,0xb4,0xc6,0xe8,0xdd,0x74,0x1f,0x4b,0xbd,0x8b,0x8a,
    0x70,0x3e,0xb5,0x66,0x48,0x03,0xf6,0x0e,0x61,0x35,0x57,0xb9,0x86,0xc1,0x1d,0x9e,
    0xe1,0xf8,0x98,0x11,0x69,0xd9,0x8e,0x94,0x9b,0x1e,0x87,0xe9,0xce,0x55,0x28,0xdf,
    0x8c,0xa1,0x89,0x0d,0xbf,0xe6,0x42,0x68,0x41,0x99,0x2d,0x0f,0xb0,0x54,0xbb,0x16,
]

INV_S_BOX = [0] * 256
for i, v in enumerate(S_BOX):
    INV_S_BOX[v] = i

# Round constants used in key expansion (Rcon)
RCON = [
    0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,0x1B,0x36,0x6C,0xD8,0xAB,0x4D
]

Nb = 4      # block size in 32-bit words (always 4 for AES)
Nk = 4      # key length in 32-bit words (4 = AES-128)
Nr = 10     # number of rounds (10 = AES-128)


# ---------------------------------------------------------------------------
# GF(2^8) multiplication used by MixColumns / InvMixColumns
# ---------------------------------------------------------------------------
def gmul(a, b):
    """Multiply two bytes in GF(2^8) with the AES reduction polynomial x^8+x^4+x^3+x+1 (0x11B)."""
    p = 0
    for _ in range(8):
        if b & 1:
            p ^= a
        hi_bit_set = a & 0x80
        a = (a << 1) & 0xFF
        if hi_bit_set:
            a ^= 0x1B
        b >>= 1
    return p & 0xFF


# ---------------------------------------------------------------------------
# Key expansion (Rijndael key schedule): 16-byte key -> 11 round keys (176 bytes)
# ---------------------------------------------------------------------------
def key_expansion(key_bytes):
    assert len(key_bytes) == 4 * Nk
    w = [list(key_bytes[4*i:4*i+4]) for i in range(Nk)]

    for i in range(Nk, Nb * (Nr + 1)):
        temp = list(w[i - 1])
        if i % Nk == 0:
            # RotWord
            temp = temp[1:] + temp[:1]
            # SubWord
            temp = [S_BOX[b] for b in temp]
            # XOR with Rcon
            temp[0] ^= RCON[i // Nk - 1]
        w.append([w[i - Nk][j] ^ temp[j] for j in range(4)])

    # Group into 11 round keys of 16 bytes each
    round_keys = []
    for r in range(Nr + 1):
        rk = []
        for c in range(4):
            rk += w[r * 4 + c]
        round_keys.append(rk)
    return round_keys


# ---------------------------------------------------------------------------
# State helpers: AES operates on a 4x4 byte matrix (column-major) called "state"
# ---------------------------------------------------------------------------
def bytes_to_state(data):
    """16 bytes -> 4x4 state matrix, filled column by column."""
    state = [[0] * 4 for _ in range(4)]
    for i in range(16):
        state[i % 4][i // 4] = data[i]
    return state


def state_to_bytes(state):
    data = [0] * 16
    for i in range(16):
        data[i] = state[i % 4][i // 4]
    return bytes(data)


# ---------------------------------------------------------------------------
# Forward round transformations
# ---------------------------------------------------------------------------
def sub_bytes(state):
    for r in range(4):
        for c in range(4):
            state[r][c] = S_BOX[state[r][c]]


def inv_sub_bytes(state):
    for r in range(4):
        for c in range(4):
            state[r][c] = INV_S_BOX[state[r][c]]


def shift_rows(state):
    for r in range(1, 4):
        state[r] = state[r][r:] + state[r][:r]


def inv_shift_rows(state):
    for r in range(1, 4):
        state[r] = state[r][-r:] + state[r][:-r]


def mix_columns(state):
    for c in range(4):
        a0, a1, a2, a3 = (state[0][c], state[1][c], state[2][c], state[3][c])
        state[0][c] = gmul(a0, 2) ^ gmul(a1, 3) ^ a2 ^ a3
        state[1][c] = a0 ^ gmul(a1, 2) ^ gmul(a2, 3) ^ a3
        state[2][c] = a0 ^ a1 ^ gmul(a2, 2) ^ gmul(a3, 3)
        state[3][c] = gmul(a0, 3) ^ a1 ^ a2 ^ gmul(a3, 2)


def inv_mix_columns(state):
    for c in range(4):
        a0, a1, a2, a3 = (state[0][c], state[1][c], state[2][c], state[3][c])
        state[0][c] = gmul(a0, 14) ^ gmul(a1, 11) ^ gmul(a2, 13) ^ gmul(a3, 9)
        state[1][c] = gmul(a0, 9) ^ gmul(a1, 14) ^ gmul(a2, 11) ^ gmul(a3, 13)
        state[2][c] = gmul(a0, 13) ^ gmul(a1, 9) ^ gmul(a2, 14) ^ gmul(a3, 11)
        state[3][c] = gmul(a0, 11) ^ gmul(a1, 13) ^ gmul(a2, 9) ^ gmul(a3, 14)


def add_round_key(state, round_key):
    for c in range(4):
        for r in range(4):
            state[r][c] ^= round_key[4 * c + r]


# ---------------------------------------------------------------------------
# Single-block encryption / decryption (AES-128, one 16-byte block)
# ---------------------------------------------------------------------------
def aes_encrypt_block(plaintext_block, key_bytes, verbose=False):
    assert len(plaintext_block) == 16
    round_keys = key_expansion(key_bytes)
    state = bytes_to_state(plaintext_block)

    add_round_key(state, round_keys[0])
    if verbose:
        print(f"Round  0 (AddRoundKey only): {state_to_bytes(state).hex().upper()}")

    for rnd in range(1, Nr):
        sub_bytes(state)
        shift_rows(state)
        mix_columns(state)
        add_round_key(state, round_keys[rnd])
        if verbose:
            print(f"Round {rnd:2d}: {state_to_bytes(state).hex().upper()}")

    # Final round (no MixColumns)
    sub_bytes(state)
    shift_rows(state)
    add_round_key(state, round_keys[Nr])
    if verbose:
        print(f"Round {Nr:2d} (final, no MixColumns): {state_to_bytes(state).hex().upper()}")

    return state_to_bytes(state)


def aes_decrypt_block(ciphertext_block, key_bytes, verbose=False):
    assert len(ciphertext_block) == 16
    round_keys = key_expansion(key_bytes)
    state = bytes_to_state(ciphertext_block)

    add_round_key(state, round_keys[Nr])
    if verbose:
        print(f"Round  0 (AddRoundKey only): {state_to_bytes(state).hex().upper()}")

    for rnd in range(Nr - 1, 0, -1):
        inv_shift_rows(state)
        inv_sub_bytes(state)
        add_round_key(state, round_keys[rnd])
        inv_mix_columns(state)
        if verbose:
            print(f"Round {Nr - rnd:2d}: {state_to_bytes(state).hex().upper()}")

    # Final round (no InvMixColumns)
    inv_shift_rows(state)
    inv_sub_bytes(state)
    add_round_key(state, round_keys[0])
    if verbose:
        print(f"Round {Nr:2d} (final, no InvMixColumns): {state_to_bytes(state).hex().upper()}")

    return state_to_bytes(state)


# ---------------------------------------------------------------------------
# PKCS#7 padding helpers (needed for arbitrary-length plaintext in ECB mode)
# ---------------------------------------------------------------------------
def pkcs7_pad(data, block_size=16):
    pad_len = block_size - (len(data) % block_size)
    return data + bytes([pad_len]) * pad_len


def pkcs7_unpad(data):
    pad_len = data[-1]
    if pad_len < 1 or pad_len > 16:
        raise ValueError("Invalid padding")
    if data[-pad_len:] != bytes([pad_len]) * pad_len:
        raise ValueError("Invalid padding")
    return data[:-pad_len]


# ---------------------------------------------------------------------------
# ECB-mode encrypt/decrypt for arbitrary-length messages (multi-block)
# ---------------------------------------------------------------------------
def aes_encrypt_ecb(plaintext_bytes, key_bytes):
    padded = pkcs7_pad(plaintext_bytes)
    out = b""
    for i in range(0, len(padded), 16):
        out += aes_encrypt_block(padded[i:i+16], key_bytes)
    return out


def aes_decrypt_ecb(ciphertext_bytes, key_bytes):
    out = b""
    for i in range(0, len(ciphertext_bytes), 16):
        out += aes_decrypt_block(ciphertext_bytes[i:i+16], key_bytes)
    return pkcs7_unpad(out)


# ---------------------------------------------------------------------------
# Demo (single test)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # --- Plaintext and key setup ---
    plaintext_text = "ovisarkar16chars"          # 16-character plaintext (128 bits)
    key_hex = "ABC231031402EFABABC231031402EFAB"[:32]  # 32-hex-digit key (128 bits)

    assert len(plaintext_text) == 16, "Plaintext must be exactly 16 characters"
    assert len(key_hex) == 32, "Key must be exactly 32 hex digits"

    plaintext_bytes = plaintext_text.encode("utf-8")
    plaintext_hex = plaintext_bytes.hex()
    key_bytes = bytes.fromhex(key_hex)

    print("Plaintext (text):", plaintext_text)
    print("Plaintext (hex): ", plaintext_hex.upper())
    print("Key (hex):       ", key_hex.upper())

    # --- Key schedule: generate and display round keys 0 to 10 ---
    round_keys = key_expansion(key_bytes)
    print("\n--- Round Keys (Key Schedule) ---")
    for i, rk in enumerate(round_keys):
        print(f"Round Key {i:2d}: {bytes(rk).hex().upper()}")

    # --- Encryption ---
    print("\n--- Encryption ---")
    ciphertext = aes_encrypt_block(plaintext_bytes, key_bytes, verbose=True)
    print("\nCiphertext (hex):", ciphertext.hex().upper())

    # --- Decryption ---
    print("\n--- Decryption ---")
    recovered_bytes = aes_decrypt_block(ciphertext, key_bytes, verbose=True)
    recovered_text = recovered_bytes.decode("utf-8")
    print("\nRecovered plaintext (hex): ", recovered_bytes.hex().upper())
    print("Recovered plaintext (text):", recovered_text)
    print("MATCH" if recovered_text == plaintext_text else "MISMATCH")