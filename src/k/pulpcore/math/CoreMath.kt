package k.pulpcore.math

import java.util.*

/**

    The CoreMath class contains fixed-point arithmetic functions and other
    useful math functions.
    <p>
    Fixed-point numbers can be used in place of floating-point numbers.
    Regarding fixed-point numbers:
    <p>
    Addition and subtraction of two fixed-point numbers is done normally, using
    + and - operators.
    <p>
    Multiplying a fixed-point number by an integer is done normally,
    using the * operator. The result is a fixed-point number.
    <p>
    Dividing a fixed-point number by an integer (the fixed-point number is the
    numerator, and the integer is the denominator) is done normally, using the /
    operator. The result is a fixed-point number.

 */
object CoreMath {

    @JvmField
    val FRACTION_BITS = 16

    @JvmField
    val FRACTION_MASK = (1 shl FRACTION_BITS) - 1

    @JvmField
    val ONE = 1 shl FRACTION_BITS

    @JvmField
    val ONE_HALF = ONE shr 1

    @JvmField
    val PI = Math.round(Math.PI * ONE.toDouble()).toInt()

    @JvmField
    val TWO_PI = Math.round(2.0 * Math.PI * ONE.toDouble()).toInt()

    @JvmField
    val ONE_HALF_PI = Math.round(0.5 * Math.PI * ONE.toDouble()).toInt()

    @JvmField
    val E = Math.round(Math.E * ONE.toDouble()).toInt()

    @JvmField
    val MAX_VALUE = Int.MAX_VALUE // (1 shl 31) -1

    @JvmField
    val MIN_VALUE = Int.MIN_VALUE // -(1 shl 31)

    /** The maximum integer value that a 32-bit fixed-point value can represent.  */
    @JvmField
    val MAX_INT_VALUE = Short.MAX_VALUE.toInt() // (1 shl 31 - FRACTION_BITS) - 1

    /** The minimum integer value that a 32-bit fixed-point value can represent.  */
    @JvmField
    val MIN_INT_VALUE = Short.MIN_VALUE.toInt() // -(1 shl 31 - FRACTION_BITS)

    /** The maximum 32-bit floating-point value that a 32-bit fixed-point value can represent.  */
    // Math.round(MAX_FLOAT_VALUE * ONE) = MAX_VALUE
    @JvmField
    val MAX_FLOAT_VALUE = Short.MAX_VALUE.toFloat() // 32768f

    /** The minimum 32-bit floating-point value that a 32-bit fixed-point value can represent.  */
    @JvmField
    val MIN_FLOAT_VALUE = -Short.MAX_VALUE.toFloat() // -32768f

    /** The maximum 64-bit floating-point value that a 32-bit fixed-point value can represent.  */
    // Found this by trial and error
    // Math.round(MAX_DOUBLE_VALUE * ONE) = MAX_VALUE
    @JvmField
    val MAX_DOUBLE_VALUE = 32767.999992370602

    /** The minimum 64-bit floating-point value that a 32-bit fixed-point value can represent.  */
    @JvmField
    val MIN_DOUBLE_VALUE = -Short.MAX_VALUE.toDouble() // -32768.0

    // For more accurate results for sine/cosine. This number was found by trial and error.
    private val TWO_PI_ERROR_FACTOR = 6

    //private static final int TWO_PI_ERROR = -11;
    private val TWO_PI_ERROR = Math.round((2.0 * Math.PI * (1 shl TWO_PI_ERROR_FACTOR).toDouble() -
            ((TWO_PI shl TWO_PI_ERROR_FACTOR) / ONE).toDouble()) * ONE).toInt()

    /** Number of fractional bits used for some internal calculations  */
    private val INTERNAL_BITS = 24

    private val rand = Random()

    /**
     * Converts an integer to a fixed-point value.
     */
    fun toFixed(n: Int) = when {
        n > MAX_INT_VALUE -> MAX_VALUE
        n < MIN_INT_VALUE -> MIN_VALUE
        else -> n shl FRACTION_BITS
    }

    /**
     * Converts a float to a fixed-point value.
     */
    fun toFixed(n: Float) = when {
        n > MAX_FLOAT_VALUE -> MAX_VALUE
        n < MIN_FLOAT_VALUE -> MIN_VALUE
        else -> Math.round(n * ONE.toFloat())
    }

    /**
     * Converts a double-precision float to a fixed-point value.
     */
    fun toFixed(n: Double) = when {
        n > MAX_DOUBLE_VALUE -> MAX_VALUE
        n < MIN_DOUBLE_VALUE -> MIN_VALUE
        else -> Math.round(n * ONE.toDouble()).toInt()
    }

    /**
     * Converts a String representing a double-precision float into a fixed-point value.
     * @throws NumberFormatException if the string is not a valid representation of a number.
     */
    fun toFixed(n: String): Int = toFixed(n.toDouble())

    /**
     * Converts a fixed-point value to a float.
     */
    fun toFloat(f: Int) = f.toFloat() / ONE

    /**
     * Converts a fixed-point value to a double.
     */
    fun toDouble(f: Int) = f.toDouble() / ONE

    /**
     * Converts a fixed-point value to an integer. Same behavior as casting a float to an int.
     */
    fun toInt(f: Int) = when {
        f < 0 -> toIntCeil(f)
        else -> toIntFloor(f)
    }

    /**
     * Converts a fixed-point value to an integer.
     */
    fun toIntFloor(f: Int) = f shr FRACTION_BITS

    /**
     * Converts a fixed-point value to an integer.
     */
    fun toIntRound(f: Int) = toIntFloor(f + ONE_HALF)

    /**
     * Converts a fixed-point value to an integer.
     */
    fun toIntCeil(f: Int) = -toIntFloor(-f)

    /**
     * Returns the fractional part of a fixed-point value (removes the integer part).
     */
    fun fracPart(f: Int) = abs(f) and FRACTION_MASK

    /**
     * Returns the floor of a fixed-point value (removes the fractional part).
     */
    fun floor(f: Int) = f and FRACTION_MASK.inv()

    /**
     * Returns the ceil of a fixed-point value.
     */
    fun ceil(f: Int) = -floor(-f)

    /**
     * Returns the fixed-point value rounded to the nearest integer location.
     */
    fun round(f: Int) = floor(f + ONE_HALF)

    /**
     * Converts a fixed-point number to a base-10 string representation.
     */
    fun toString(f: Int) = formatNumber(
            abs(toInt(f)),
            fracPart(f) shl 32 - FRACTION_BITS, f < 0,
            1, 7, false)

    /**
     * Converts a fixed-point number to a base-10 string representation using
     * the specified number of fractional digits.
     */
    fun toString(f: Int, numFractionalDigits: Int) = formatNumber(
            abs(toInt(f)),
            fracPart(f) shl 32 - FRACTION_BITS, f < 0,
            numFractionalDigits, numFractionalDigits, false)

    /**
     * Converts a fixed-point number to a base-10 string representation.
     * @param f the fixed-point number
     * @param minFracDigits the minimum number of digits to show after
     * the decimal point.
     * @param maxFracDigits the maximum number of digits to show after
     * the decimal point.
     * @param grouping if (true, uses the grouping character (',')
     * to seperate groups in the integer portion of the number.
     */
    fun toString(f: Int, minFracDigits: Int, maxFracDigits: Int, grouping: Boolean) = formatNumber(
            abs(toInt(f)),
            fracPart(f) shl 32 - FRACTION_BITS, f < 0,
            minFracDigits, maxFracDigits, grouping)

    /**
     * Converts an integer to a base-10 string representation.
     */
    fun intToString(n: Int) = formatNumber(abs(n), 0, n < 0, 0, 0, false)

    /**
     * Converts a integer to a base-10 string representation using
     * the specified number of fractional digits.
     */
    fun intToString(n: Int, numFractionalDigits: Int) = formatNumber(
            abs(n),
            0, n < 0,
            numFractionalDigits, numFractionalDigits, false)

    /**
     * Converts an integer to a base-10 string representation.
     * @param n the integer
     * @param minFracDigits the minimum number of digits to show after
     * the decimal point.
     * @param maxFracDigits the maximum number of digits to show after
     * the decimal point.
     * @param grouping if (true, uses the grouping character (',')
     * to seperate groups in the integer portion of the number.
     */
    fun intToString(n: Int, minFracDigits: Int, maxFracDigits: Int, grouping: Boolean) = formatNumber(
            abs(n),
            0, n < 0,
            minFracDigits, maxFracDigits, grouping)

    // TODO review all 2, e.g. intPart2 --> intPart
    /**
     * Converts a number to a base-10 string representation.
     * @param intPart2 the integer part of the number
     * @param fracPart the fractional part, a 32-bit fixed point value.
     * @param minFracDigits the minimum number of digits to show after
     * the decimal point.
     * @param maxFracDigits the maximum number of digits to show after
     * the decimal point.
     * @param intPartGrouping if (true, uses the groupong character (',')
     * to seperate groups in the integer portion of the number.
     */
    private fun formatNumber(intPart2: Int, fracPart: Int, negative: Boolean,
                             minFracDigits: Int, maxFracDigits: Int, intPartGrouping: Boolean): String {
        var intPart = intPart2
        val buffer = StringBuffer()
        val one = 1L shl 32
        val mask = one - 1
        var frac = fracPart.toLong() and mask

        // Round up if needed
        if (maxFracDigits < 10) {
            var place: Long = 1
            for (i in 0 until maxFracDigits) {
                place *= 10
            }
            frac += (1L shl 31) / place
            if (frac >= one) {
                intPart++
            }
        }

        // Convert integer part
        if (!intPartGrouping || intPart == 0) {
            if (negative) {
                buffer.append('-')
            }
            buffer.append(intPart)
        } else {
            var i = 0
            while (intPart > 0) {
                if (i == 3) {
                    buffer.insert(0, ',')
                    i = 0
                }
                val ch = (intPart % 10 + '0'.toInt()).toChar()
                buffer.insert(0, ch)
                intPart /= 10
                i++
            }

            if (negative) {
                buffer.insert(0, '-')
            }
        }

        if (maxFracDigits == 0 || fracPart == 0 && minFracDigits == 0) {
            return buffer.toString()
        }

        buffer.append('.')

        // Convert fractional part
        var numFracDigits = 0
        while (true) {
            frac = (frac and mask) * 10L

            if (frac == 0L) {
                buffer.append('0')
            } else {
                buffer.append(('0'.toLong() + frac.ushr(32) % 10).toChar())
            }

            numFracDigits++
            if (numFracDigits == maxFracDigits || frac == 0L && numFracDigits >= minFracDigits) {
                break
            }
        }

        // Remove uneeded trailing zeros
        if (numFracDigits > minFracDigits) {
            val len = numFracDigits - minFracDigits
            for (i in 0..len - 1) {
                if (buffer[buffer.length - 1] == '0') {
                    buffer.setLength(buffer.length - 1)
                } else {
                    break
                }
            }
        }

        return buffer.toString()
    }

    //
    // Bit manipulation
    //

    /**
     * Returns true if the number (greater than 1) is a power of two.
     */
    fun isPowerOfTwo(n: Int) = n and n - 1 == 0

    /**
     * Counts the number of "on" bits in an integer.
     */
    fun countBits(n: Int): Int {
        /*
        int count = 0;
        while (n > 0) {
            count += (n & 1);
            n >>= 1;
        }
        return count;
        */
        var count = n

        count = (count shr 1  and 0x55555555) + (count and 0x55555555)
        count = (count shr 2  and 0x33333333) + (count and 0x33333333)
        count = (count shr 4  and 0x0F0F0F0F) + (count and 0x0F0F0F0F)
        count = (count shr 8  and 0x00FF00FF) + (count and 0x00FF00FF)
        count = (count shr 16 and 0x0000FFFF) + (count and 0x0000FFFF)

        return count
    }

    /**
     * Returns the log base 2 of an integer greater than 0. The returned value
     * is equal to `Math.floor(Math.log(n) / Math.log(2))`.
     */
    fun log2(n2: Int): Int {
        var n = n2
        //if (n <= 1) {
        //    throw new ArithmeticException("NaN");
        //}
        var count = 0
        while (true) {
            n = n shr 1
            if (n == 0) {
                return count
            }
            count++
        }

        /*
        int count = 0;

        if ((n & 0xFFFF0000) != 0) {
            n >>= 16;
            count = 16;
        }
        if ((n & 0xFF00) != 0) {
            n >>= 8;
            count |= 8;
        }
        if ((n & 0xF0) != 0) {
            n >>= 4;
            count |= 4;
        }
        if ((n & 0xC) != 0) {
            n >>= 2;
            count |= 2;
        }
        if ((n & 0x2) != 0) {
            //n >>= 1;
            count |= 1;
        }

        return count;
*/
    }

    //
    // Integer math
    //

    /**
     * Clamps a number between two values. If the number <= min returns min; if the number >= max
     * returns max; otherwise returns the number.
     */
    fun clamp(n: Int, min: Int, max: Int) = when {
        n <= min -> min
        n >= max -> max
        else -> n
    }

    /**
     * Clamps a number between two values. If the number <= min returns min; if the number >= max
     * returns max; otherwise returns the number.
     */
    fun clamp(n: Float, min: Float, max: Float) = when {
        n <= min -> min
        n >= max -> max
        else -> n
    }

    /**
     * Clamps a number between two values. If the number <= min returns min; if the number >= max
     * returns max; otherwise returns the number.
     */
    fun clamp(n: Double, min: Double, max: Double) = when {
        n <= min -> min
        n >= max -> max
        else -> n
    }

    /**
     * Returns the sign of a number.
     */
    fun sign(n: Int) = when {
        n > 0 -> 1
        n < 0 -> -1
        else -> 0
    }

    /**
     * Returns the sign of a number.
     */
    fun sign(n: Double) = when {
        n > 0 -> 1
        n < 0 -> -1
        else -> 0
    }

    /**
     * Returns the absolute value of a number.
     */
    fun abs(n: Int) = when {
        n >= 0 -> n
        else -> -n
    }

    /**
     * Divides the number, n, by the divisor, d, rounding the result to the
     * nearest integer.
     */
    fun intDivRound(n: Int, d: Int) = when {
        ((d > 0) xor (n > 0)) -> (n - (d shr 1)) / d
        else -> (n + (d shr 1)) / d
    }

    /**
     * Divides the number, n, by the divisor, d, returning the nearest integer
     * less than or equal to the result.
     */
    fun intDivFloor(n: Int, d: Int) = when {
        d > 0 -> when {
            n < 0 -> (n - d + 1) / d
            else -> n / d
        }
        d < 0 -> when {
            n > 0 -> (n - d - 1) / d
            else -> n / d
        }
        else -> n / d
    }

    /**
     * Divides the number, n, by the divisor, d, returning the nearest integer
     * greater than or equal to the result.
     */
    fun intDivCeil(n: Int, d: Int) = -intDivFloor(-n, d)

    //
    // Fixed-point math
    //

    /**
     * Multiplies two fixed-point numbers together.
     */
    fun mul(f1: Int, f2: Int) = (f1.toLong() * f2 shr FRACTION_BITS).toInt()

    /**
     * Multiplies two fixed-point numbers together.
     */
    fun mul(f1: Long, f2: Long) = f1 * f2 shr FRACTION_BITS

    /**
     * Divides the first fixed-point number by the second fixed-point number.
     */
    fun div(f1: Int, f2: Int) = ((f1.toLong() shl FRACTION_BITS) / f2).toInt()

    /**
     * Divides the first fixed-point number by the second fixed-point number.
     */
    fun div(f1: Long, f2: Long) = (f1 shl FRACTION_BITS) / f2

    /**
     * Multiplies the first two fixed-point numbers together, then divides by
     * the third fixed-point number.
     */
    fun mulDiv(f1: Int, f2: Int, f3: Int) = (f1.toLong() * f2 / f3).toInt()

    /**
     * Multiplies the first two fixed-point numbers together, then divides by
     * the third fixed-point number.
     */
    fun mulDiv(f1: Long, f2: Long, f3: Long) = f1 * f2 / f3

    //
    // Logs and powers
    //

    fun sqrt(fx2: Int): Int {
        var fx = fx2
        if (fx < 0) {
            throw ArithmeticException("NaN")
        }

        if (fx == 0 || fx == ONE) {
            return fx
        }

        // invert numbers less than one (if they aren't too small)
        var invert = false
        if (fx < ONE && fx > 6) {
            invert = true
            fx = div(ONE, fx)
        }

        var iterations = 16
        if (fx > ONE) {
            // number of iterations == (number of bits in number) / 2
            var s = fx
            iterations = 0
            while (s > 0) {
                s = s shr 2
                iterations++
            }
        }

        // Newton's iteration
        var l = (fx shr 1) + 1
        for (i in 1 until iterations) {
            l = l + div(fx, l) shr 1
        }

        // undo the inversion
        return if (invert) {
            div(ONE, l)
        } else l

    }

    fun sqrt(fx2: Long): Long {
        var fx = fx2
        if (fx < 0L) {
            throw ArithmeticException("NaN")
        }

        if (fx == 0L || fx == ONE.toLong()) {
            return fx
        }

        // Invert numbers less than one (if they aren't too small)
        var invert = false
        if (fx < ONE && fx > 6L) {
            invert = true
            fx = div(ONE.toLong(), fx)
        }

        var iterations = 16
        if (fx > ONE) {
            // Number of iterations == (number of bits in number) / 2
            var s = fx
            iterations = 0
            while (s > 0L) {
                s = s shr 2
                iterations++
            }
        }

        // Newton's iteration
        var l = (fx shr 1) + 1
        for (i in 1..iterations - 1) {
            l = l + div(fx, l) shr 1
        }

        // Undo the inversion
        return if (invert) {
            div(ONE.toLong(), l)
        } else l

    }

    fun dist(x1: Int, y1: Int, x2: Int, y2: Int): Long {
        val dx = (x1 - x2).toLong()
        val dy = (y1 - y2).toLong()
        return sqrt(mul(dx, dx) + mul(dy, dy))
    }

    //
    // Fixed-point Trigonometry
    //

    // TODO fx2 -> fx
    /**
     * Returns the sine of the specified fixed-point radian value.
     */
    fun sin(fx2: Int): Int {
        var fx = fx2

        if (fx == 0) {
            return 0
        }

        // reduce range to -2*pi and 2*pi
        val s = fx / TWO_PI
        if (abs(s) >= 1 shl TWO_PI_ERROR_FACTOR) {
            // fix any error for large values of fx
            fx -= s * TWO_PI + (s shr TWO_PI_ERROR_FACTOR) * TWO_PI_ERROR
        } else {
            fx -= s * TWO_PI
        }

        // reduce range to -pi/2 and pi/2
        // this allows us to limit the number of iterations in the maclaurin series
        if (fx > PI) {
            fx -= -TWO_PI
        } else if (fx < -PI) {
            fx += TWO_PI
        }
        if (fx > ONE_HALF_PI) {
            fx = PI - fx
        } else if (fx < -ONE_HALF_PI) {
            fx = -PI - fx
        }

        // Helps with rotation appearance near 90, 180, 270, 360, etc.
        if (abs(fx) < 32) {
            return 0
        } else if (abs(fx - ONE_HALF_PI) < 32) {
            return ONE
        } else if (abs(fx + ONE_HALF_PI) < 32) {
            return -ONE
        }

        // Maclaurin power series
        val fxSquared = mul(fx, fx)
        val d = mul((1 shl INTERNAL_BITS) / (2 * 3 * 4 * 5 * 6 * 7 * 8 * 9), fxSquared)
        val c = mul(d - (1 shl INTERNAL_BITS) / (2 * 3 * 4 * 5 * 6 * 7), fxSquared)
        val b = mul(c + (1 shl INTERNAL_BITS) / (2 * 3 * 4 * 5), fxSquared)
        val a = mul(b - (1 shl INTERNAL_BITS) / (2 * 3), fxSquared)
        val sine = mul(a + (1 shl INTERNAL_BITS), fx)
        return sine shr INTERNAL_BITS - FRACTION_BITS
    }

    /**
     * Returns the cosine of the specified fixed-point radian value.
     */
    fun cos(fx: Int) = when {
        fx == 0 -> ONE
        fx < 0 -> sin(ONE_HALF_PI - TWO_PI - fx) // make up for potential overflow
        else -> sin(ONE_HALF_PI - fx)
    }

    /**
     * Returns the tangent of the specified fixed-point radian value.
     */
    fun tan(fx: Int): Int {
        val cos = cos(fx)
        when {
            cos == 0 -> MAX_VALUE
        }
        return if (cos == 0) MAX_VALUE else div(sin(fx), cos)
    }

    /**
     * Returns the cotangent of the specified fixed-point radian value.
     */
    fun cot(fx: Int): Int {
        val sin = sin(fx)
        return if (sin == 0) {
            MAX_VALUE
        } else {
            div(cos(fx), sin)
        }
    }

    /**
     * Returns the arcsine of the specified fixed-point value.
     */
    fun asin(fx: Int): Int {
        return if (abs(fx) > ONE) {
            throw ArithmeticException("NaN")
        } else if (fx == ONE) {
            ONE_HALF_PI
        } else if (fx == -ONE) {
            -ONE_HALF_PI
        } else {
            atan(div(fx, sqrt(ONE - mul(fx, fx))))
        }
    }

    /**
     * Returns the arccosine of the specified fixed-point value.
     */
    fun acos(fx: Int) = ONE_HALF_PI - asin(fx)

    /**
     * Returns the arctangent of the specified fixed-point value.
     */
    fun atan(fx2: Int): Int {
        var fx = fx2
        var negative = false
        var invert = false
        if (fx == 0) {
            return 0
        }
        if (fx < 0) {
            negative = true
            fx = -fx
        }

        // Avoid overflow
        if (fx > ONE) {
            invert = true
            fx = div(ONE, fx)
        }

        // Approximation from Ranko at http://www.lightsoft.co.uk/PD/stu/stuchat37.html
        // r(x) = (x + 0.43157974*x^3)/(1 + 0.76443945*x^2 + 0.05831938*x^4)
        val fxPow2 = mul(fx, fx)
        val fxPow3 = mul(fxPow2, fx)
        val fxPow4 = mul(fxPow3, fx)
        val numer = fx + mul(28284, fxPow3)
        val denom = ONE + mul(50098, fxPow2) + mul(3822, fxPow4)
        var answer = div(numer, denom)

        if (invert) {
            answer = ONE_HALF_PI - answer
        }
        if (negative) {
            answer = -answer
        }
        return answer
    }

    /**
     * Returns in the range from -pi to pi.
     */
    fun atan2(fy: Int, fx: Int): Int {
        if (fy == 0) {
            return if (fx < 0) {
                PI
            } else {
                0
            }
        } else if (fx == 0) {
            return if (fy < 0) {
                -ONE_HALF_PI
            } else {
                ONE_HALF_PI
            }
        } else {
            val answer = atan(abs(div(fy, fx)))
            return if (fy > 0 && fx < 0) {
                PI - answer
            } else if (fy < 0 && fx < 0) {
                answer - PI
            } else if (fy < 0 && fx > 0) {
                -answer
            } else {
                answer
            }
        }
    }

    //
    // Random number generation and Noise functions
    //

    /**
     * Returns a random integer from 0 to max, inclusive
     */
    fun rand(max: Int): Int {
        return rand(0, max)
    }

    /**
     * Returns a random integer from min to max, inclusive
     */
    fun rand(min: Int, max: Int): Int {
        val range = max.toLong() - min.toLong() + 1
        val r = rand.nextInt().toLong() and 0xffffffffL
        val value = (min + r * range / 0xffffffffL).toInt()

        // Bounds check is probably not needed.
        return when {
            value < min -> min
            value > max -> max
            else -> value
        }
    }

    /**
     * Returns a random double from 0 to max, inclusive
     */
    fun rand(max: Double) = rand(0.0, max)

    /**
     * Returns a random double from min to max, inclusive
     */
    fun rand(min: Double, max: Double): Double {
        val value = min + rand.nextDouble() * (max - min)
        // Bounds check is probably not needed.
        return when {
            value < min -> min
            value > max -> max
            else -> value
        }
    }

    /**
     * Returns a random boolean.
     */
    fun rand() = rand(0, 1) == 0

    /**
     * Returns true if a random event occurs.
     * @param percent The probability of the event occuring, from 0 (never) to 100 (always).
     */
    fun randChance(percent: Int) = rand(1, 100) <= percent

    /**
     * @return a 32-bit signed integer
     */
    fun noise(n2: Int): Int {
        var n = n2

        // A common noise function.
        // Note: these numbers are all primes.
        n = n shl 13 xor n
        return n * (n * n * 15731 + 789221) + 1376312589
    }

    /**
     * @return a 32-bit signed integer
     */
    fun noise(x: Int, y: Int) = noise(x + y * 57)

    /**
     * @return a 32-bit signed integer
     */
    fun smoothNoise(x: Int) = ((noise(x).toLong() shr 1) + (noise(x - 1).toLong() shr 2) + (noise(x + 1).toLong() shr 2)).toInt()

    /**
     * @return a 32-bit signed integer
     */
    fun smoothNoise(x: Int, y: Int): Int {
        val corners = noise(x - 1, y - 1).toLong() +
                noise(x + 1, y - 1).toLong() +
                noise(x - 1, y + 1).toLong() +
                noise(x + 1, y + 1).toLong() shr 4
        val sides = noise(x - 1, y).toLong() +
                noise(x + 1, y).toLong() +
                noise(x, y - 1).toLong() +
                noise(x, y + 1).toLong() shr 3
        val center = (noise(x, y) shr 2).toLong()

        return (corners + sides + center).toInt()
    }

    /**
     * @return a 32-bit signed integer
     */
    fun interpolatedNoise(fx: Int): Int {
        val x = fx shr FRACTION_BITS
        val f = fx and FRACTION_MASK

        if (f == 0) {
            return smoothNoise(x)
        } else {
            val n1 = smoothNoise(x)
            val n2 = smoothNoise(x + 1)

            return cosineInterpolate(n1, n2, f)
        }
    }

    /**
     * @return a 32-bit signed integer
     */
    fun interpolatedNoise(fx: Int, fy: Int): Int {
        val x = fx shr FRACTION_BITS
        val y = fy shr FRACTION_BITS

        val v1 = smoothNoise(x, y)
        val v2 = smoothNoise(x + 1, y)
        val v3 = smoothNoise(x, y + 1)
        val v4 = smoothNoise(x + 1, y + 1)

        val n1 = cosineInterpolate(v1, v2, fx and FRACTION_MASK)
        val n2 = cosineInterpolate(v3, v4, fx and FRACTION_MASK)

        return cosineInterpolate(n1, n2, fy and FRACTION_MASK)
    }

    /**
     * @param fx fixed-point value
     * @param persistence fixed-point value <= 1.
     * @return a 32-bit signed integer
     *
     * From http://freespace.virgin.net/hugo.elias/models/m_perlin.htm
     */
    fun perlinNoise(fx: Int, persistence: Int, numOctaves: Int): Int {
        var total: Long = 0

        var amplitude = ONE - persistence

        for (i in 0..numOctaves - 1) {

            total += mul(interpolatedNoise(fx shl i), amplitude).toLong()

            amplitude = mul(amplitude, persistence)

        }

        return total.toInt()
    }

    /**
     * @param fx fixed-point value
     * @param fy fixed-point value
     * @param persistence fixed-point value <= 1.
     * @return a 32-bit signed integer
     *
     * From http://freespace.virgin.net/hugo.elias/models/m_perlin.htm
     */
    fun perlinNoise(fx: Int, fy: Int, persistence: Int, numOctaves: Int): Int {
        var total: Long = 0

        var amplitude = ONE - persistence

        for (i in 0..numOctaves - 1) {

            total += mul(interpolatedNoise(fx shl i, fy shl i), amplitude).toLong()

            amplitude = mul(amplitude, persistence)

        }

        return total.toInt()
    }

    //
    // Interpolation
    //

    /**
     * Performs a 1-dimensional linear interpolation between values n1 and n2. The f
     * parameter is a fixed-point value from 0.0 to 1.0.
     * If f is less than 0 or greater than 1, extrapolation is calculated.
     */
    fun interpolate(n1: Int, n2: Int, f: Int) = mul(n1, ONE - f) + mul(n2, f) // This could cause overflow errors: return n1 + mul(n2 - n1, f);

    /**
     * Performs a 1-dimensional cosine interpolation between values n1 and n2. The f
     * parameter is a fixed-point value from 0.0 to 1.0.
     * If f is less than 0 or greater than 1, extrapolation is calculated.
     */
    fun cosineInterpolate(n1: Int, n2: Int, f2: Int): Int {
        var f = f2

        f = mul(f, PI)
        f = ONE - cos(f) shr 1

        return interpolate(n1, n2, f)
    }

    // This supposedly resembles cosine interpolation by using the graph
    // y = 3x^2 - 2x^3
    fun quickCurveInterpolate(n1: Int, n2: Int, f: Int): Int {
        val fSquared = mul(f, f)
        val fCubed = mul(fSquared, f)

        val fi = ONE - f
        val fiSquared = mul(fi, fi)
        val fiCubed = mul(fiSquared, fi)

        return mul(n1, 3 * fiSquared - 2 * fiCubed) + mul(n2, 3 * fSquared - 2 * fCubed)
    }

    /**
     * Performs a 1-dimensional cubic interpolation between values n1 and n2. The f
     * parameter is a fixed-point value from 0.0 to 1.0.
     * If f is less than 0 or greater than 1, extrapolation is calculated.
     */
    fun cubicInterpolate(n0: Int, n1: Int, n2: Int, n3: Int, f: Int): Int {
        val f2 = mul(f, f)
        val f3 = mul(f2, f)

        val p = n3.toLong() - n2 - (n0.toLong() - n1)
        val q = n0.toLong() - n1 - p
        val r = n2.toLong() - n0
        val s = n1.toLong()

        return (mul(p, f3.toLong()) + mul(q, f2.toLong()) + mul(r, f.toLong()) + s).toInt()
    }

}
