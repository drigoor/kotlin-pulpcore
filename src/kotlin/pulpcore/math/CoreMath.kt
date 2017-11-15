package kotlin.pulpcore.math

const val FRACTION_BITS = 16

//const val MAX_VALUE = (1 shl 31) -1
//
//const val MIN_VALUE = -(1 shl 31)

const val MAX_INT_VALUE = (1 shl 31 - FRACTION_BITS) - 1

const val MIN_INT_VALUE = -(1 shl 31 - FRACTION_BITS)

class CoreMath private constructor() {

    /**
     * Converts an integer to a fixes-point value.
     */
    fun toFixed(n: Int) = when {
        n > MAX_INT_VALUE -> Int.MAX_VALUE
        n < MIN_INT_VALUE -> Int.MIN_VALUE
        else -> n shl FRACTION_BITS
    }

}