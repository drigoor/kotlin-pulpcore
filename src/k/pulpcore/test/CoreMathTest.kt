package k.pulpcore.test

import k.pulpcore.math.CoreMath

object CoreMathTest {

    const val NUM_TESTS = 100000

    fun sinNearZero(delta: Double = 0.0009) {
        for (i in 0..NUM_TESTS) {
            val angle = CoreMath.rand(-2 * Math.PI, 2 * Math.PI)
            val expectedResult = Math.sin(angle)
            val actualResult = CoreMath.toDouble(CoreMath.sin(CoreMath.toFixed(angle)))

            when {
                Math.abs(expectedResult.minus(actualResult)) > delta -> println("WRONG test 'sinNearZero' expected result: $expectedResult VS actual result:  $actualResult")
            }

            println("expected: $expectedResult actual: $actualResult")
        }
    }
}

fun main(args: Array<String>) {
    CoreMathTest.sinNearZero()

//    println(100.toDouble() / CoreMath.ONE.toDouble())


//    println("1st -> " + (2.0 * Math.PI * (1 shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble()))
//    // 402,1238596594935
////    println("2nd -> " + (CoreMath.TWO_PI shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble())
////    // 2.63536E7
//    println("3rd -> " + ((CoreMath.TWO_PI shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble() / CoreMath.ONE.toDouble()))
//    // 402,1240234375
//    println()
//    val result: Double = ((2.0 * Math.PI * (1 shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble()) - ((CoreMath.TWO_PI shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble() / CoreMath.ONE.toDouble())) * CoreMath.ONE.toDouble()
//    println("--> " + Math.round(((2.0 * Math.PI * (1 shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble()) - ((CoreMath.TWO_PI shl CoreMath.TWO_PI_ERROR_FACTOR).toDouble() / CoreMath.ONE.toDouble())) * CoreMath.ONE.toDouble()).toInt())
}