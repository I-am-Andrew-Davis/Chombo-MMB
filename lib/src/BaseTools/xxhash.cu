/*
 * xxHash - Extremely Fast Hash algorithm
 * Header File
 * Copyright (C) 2012-2021 Yann Collet
 *
 * BSD 2-Clause License (https://www.opensource.org/licenses/bsd-license.php)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above
 *      copyright notice, this list of conditions and the following disclaimer
 *      in the documentation and/or other materials provided with the
 *      distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * You can contact the author at:
 *   - xxHash homepage: https://www.xxhash.com
 *   - xxHash source repository: https://github.com/Cyan4973/xxHash
 */

#include <cstdint>

#include "xxhash.H"

namespace xxhash
{

/* #define instead of static const, to be used as initializers */
#define XXH_PRIME32_1  0x9E3779B1U  /*!< 0b10011110001101110111100110110001 */
#define XXH_PRIME32_2  0x85EBCA77U  /*!< 0b10000101111010111100101001110111 */
#define XXH_PRIME32_3  0xC2B2AE3DU  /*!< 0b11000010101100101010111000111101 */
#define XXH_PRIME32_4  0x27D4EB2FU  /*!< 0b00100111110101001110101100101111 */
#define XXH_PRIME32_5  0x165667B1U  /*!< 0b00010110010101100110011110110001 */

#define XXH_rotl32(x,r) (((x) << (r)) | ((x) >> (32 - (r))))

// Byte swap is the same as used for google cityhash
#if defined(__CUDACC__)

__host__ __device__ inline constexpr uint32_t
bswap_32(const uint32_t a_x)
{
  return __byte_perm(a_x, 0, 0x0123);
}

#else

#include <byteswap.h>

#endif
#define XXH_swap32 bswap_32

// Portable and safe solution. Generally efficient.
HOSTDEVICE static xxh_u32
XXH_read32(const void* memPtr)
{
    xxh_u32 val;
    std::memcpy(&val, memPtr, sizeof(val));
    return val;
}

// Normal stripe processing routine.
HOSTDEVICE static xxh_u32
XXH32_round(xxh_u32 acc, xxh_u32 input)
{
    acc += input * XXH_PRIME32_2;
    acc  = XXH_rotl32(acc, 13);
    acc *= XXH_PRIME32_1;
#if (defined(__SSE4_1__) || defined(__aarch64__)) && !defined(XXH_ENABLE_AUTOVECTORIZE)
    XXH_COMPILER_GUARD(acc);
#endif
    return acc;
}

//
HOSTDEVICE XXH_FORCE_INLINE xxh_u32
XXH_readLE32(const void* ptr)
{
    return XXH_CPU_LITTLE_ENDIAN ? XXH_read32(ptr) : XXH_swap32(XXH_read32(ptr));
}

#define XXH_get32bits(p) XXH_readLE32_align(p, align)
//
HOSTDEVICE XXH_FORCE_INLINE xxh_u32
XXH_readLE32_align(const void* ptr, XXH_alignment align)
{
    if (align==XXH_unaligned) {
        return XXH_readLE32(ptr);
    } else {
        return XXH_CPU_LITTLE_ENDIAN ? *(const xxh_u32*)ptr : XXH_swap32(*(const xxh_u32*)ptr);
    }
}

// Mixes all bits to finalize the hash.
HOSTDEVICE static xxh_u32 XXH32_avalanche(xxh_u32 h32)
{
    h32 ^= h32 >> 15;
    h32 *= XXH_PRIME32_2;
    h32 ^= h32 >> 13;
    h32 *= XXH_PRIME32_3;
    h32 ^= h32 >> 16;
    return(h32);
}

// Processes the last 0-15 bytes of @p ptr.
HOSTDEVICE static xxh_u32
XXH32_finalize(xxh_u32 h32, const xxh_u8* ptr, size_t len, XXH_alignment align)
{
#define XXH_PROCESS1 do {                           \
    h32 += (*ptr++) * XXH_PRIME32_5;                \
    h32 = XXH_rotl32(h32, 11) * XXH_PRIME32_1;      \
} while (0)

#define XXH_PROCESS4 do {                           \
    h32 += XXH_get32bits(ptr) * XXH_PRIME32_3;      \
    ptr += 4;                                   \
    h32  = XXH_rotl32(h32, 17) * XXH_PRIME32_4;     \
} while (0)

    if (ptr==NULL) XXH_ASSERT(len == 0);

    /* Compact rerolled version; generally faster */
    if (!XXH32_ENDJMP) {
        len &= 15;
        while (len >= 4) {
            XXH_PROCESS4;
            len -= 4;
        }
        while (len > 0) {
            XXH_PROCESS1;
            --len;
        }
        return XXH32_avalanche(h32);
    } else {
         switch(len&15) /* or switch(bEnd - p) */ {
           case 12:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 8:       XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 4:       XXH_PROCESS4;
                         return XXH32_avalanche(h32);

           case 13:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 9:       XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 5:       XXH_PROCESS4;
                         XXH_PROCESS1;
                         return XXH32_avalanche(h32);

           case 14:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 10:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 6:       XXH_PROCESS4;
                         XXH_PROCESS1;
                         XXH_PROCESS1;
                         return XXH32_avalanche(h32);

           case 15:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 11:      XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 7:       XXH_PROCESS4;
                         XXH_FALLTHROUGH;
           case 3:       XXH_PROCESS1;
                         XXH_FALLTHROUGH;
           case 2:       XXH_PROCESS1;
                         XXH_FALLTHROUGH;
           case 1:       XXH_PROCESS1;
                         XXH_FALLTHROUGH;
           case 0:       return XXH32_avalanche(h32);
        }
        XXH_ASSERT(0);
        return h32;   /* reaching this point is deemed impossible */
    }
}

// The implementation for @ref XXH32().
HOSTDEVICE XXH_FORCE_INLINE xxh_u32
XXH32_endian_align(const xxh_u8* input, size_t len, xxh_u32 seed, XXH_alignment align)
{
    xxh_u32 h32;

    if (input==NULL) XXH_ASSERT(len == 0);

    if (len>=16) {
        const xxh_u8* const bEnd = input + len;
        const xxh_u8* const limit = bEnd - 15;
        xxh_u32 v1 = seed + XXH_PRIME32_1 + XXH_PRIME32_2;
        xxh_u32 v2 = seed + XXH_PRIME32_2;
        xxh_u32 v3 = seed + 0;
        xxh_u32 v4 = seed - XXH_PRIME32_1;

        do {
            v1 = XXH32_round(v1, XXH_get32bits(input)); input += 4;
            v2 = XXH32_round(v2, XXH_get32bits(input)); input += 4;
            v3 = XXH32_round(v3, XXH_get32bits(input)); input += 4;
            v4 = XXH32_round(v4, XXH_get32bits(input)); input += 4;
        } while (input < limit);

        h32 = XXH_rotl32(v1, 1)  + XXH_rotl32(v2, 7)
            + XXH_rotl32(v3, 12) + XXH_rotl32(v4, 18);
    } else {
        h32  = seed + XXH_PRIME32_5;
    }

    h32 += (xxh_u32)len;

    return XXH32_finalize(h32, input, len&15, align);
}

/*!
 * @brief Calculates the 32-bit hash of @p input using xxHash32.
 *
 * Speed on Core 2 Duo @ 3 GHz (single thread, SMHasher benchmark): 5.4 GB/s
 *
 * @param input The block of data to be hashed, at least @p length bytes in size.
 * @param length The length of @p input, in bytes.
 * @param seed The 32-bit seed to alter the hash's output predictably.
 *
 * @pre
 *   The memory between @p input and @p input + @p length must be valid,
 *   readable, contiguous memory. However, if @p length is `0`, @p input may be
 *   `NULL`. In C++, this also must be *TriviallyCopyable*.
 *
 * @return The calculated 32-bit hash value.
 *
 * @see
 *    XXH64(), XXH3_64bits_withSeed(), XXH3_128bits_withSeed(), XXH128():
 *    Direct equivalents for the other variants of xxHash.
 * @see
 *    XXH32_createState(), XXH32_update(), XXH32_digest(): Streaming version.
 */
HOSTDEVICE XXH_PUBLIC_API XXH32_hash_t
XXH32 (const void* input, size_t len, XXH32_hash_t seed)
{
    if (XXH_FORCE_ALIGN_CHECK) {
        if ((((size_t)input) & 3) == 0) {   /* Input is 4-bytes aligned, leverage the speed benefit */
            return XXH32_endian_align((const xxh_u8*)input, len, seed, XXH_aligned);
    }   }

    return XXH32_endian_align((const xxh_u8*)input, len, seed, XXH_unaligned);
}

}  // namespace xxhash
