import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * This JUnit test class provides a comprehensive test suite for the Ex1 assignment.
 * It evaluates correctness, robustness, edge-case handling and mathematical accuracy
 * of all polynomial utility functions, including:
 *   - Evaluation (f)
 *   - Addition, subtraction, multiplication
 *   - Derivative and normalization
 *   - Polynomial parsing and string formatting
 *   - Root finding, area computation and curve length
 *   - Construction of polynomials from points
 *   - EPS-based equality and null-input behavior
 *
 * Over 70 tests cover normal cases, boundary conditions, invalid inputs,
 * large numbers, reversed intervals, ZERO polynomial behavior, symmetry,
 * and consistency checks between mathematical definitions.
 *
 * Each test verifies one well-defined property (unit testing),
 * ensuring that any change in the implementation will immediately detect regressions.
 */


class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}
    /**
     * Tests the derivative(double[]) function on a cubic polynomial.
     * The test verifies that the derivative of [2,4,6,8] is correctly computed as [4,12,24],
     * and uses equals() to confirm polynomial equality.
     */
    @Test
    public void testDerivative() {
        double [] po = {2,4,6,8};
        double [] poD = {4,12,24};
        double [] ans = Ex1.derivative(po);
        if (!Ex1.equals(ans,poD)){
            fail();
        }
    }
    /**
     * Tests the poly(double[]) and getPolynomFromString(String) functions together.
     * The test converts a polynomial array to a string using poly(),
     * then parses it back using getPolynomFromString(),
     * and verifies that the reconstructed polynomial equals the original.
     */
    @Test
    public void testPolyToString(){
        double [] po = {2.2,-4,6,8};
        String s = Ex1.poly(po);
        double [] test = Ex1.getPolynomFromString(s);
        if(!Ex1.equals(test,po)){
            fail();
        }
    }
    /**
     * Tests that f(null, x) throws NullPointerException.
     */
    @Test
    public void test_f_null() {
        assertThrows(NullPointerException.class, () -> Ex1.f(null, 5));
    }

    /**
     * Tests f(x) behavior when polynomial contains NaN or Infinity.
     */
    @Test
    public void test_f_nan_inf() {
        double[] p = {Double.NaN, 1};
        assertTrue(Double.isNaN(Ex1.f(p, 2)));

        double[] p2 = {Double.POSITIVE_INFINITY, 1};
        assertTrue(Double.isInfinite(Ex1.f(p2, 1)));
    }

    /**
     * Tests f(x) on a constant polynomial.
     */
    @Test
    public void testFSimpleConstant() {
        double[] p = {5}; // f(x) = 5
        assertEquals(5, Ex1.f(p, 10), Ex1.EPS);
        assertEquals(5, Ex1.f(p, -3), Ex1.EPS);
        assertEquals(5, Ex1.f(p, 0), Ex1.EPS);
    }

    /**
     * Tests f(x) on a linear polynomial.
     */
    @Test
    public void testFLinear() {
        double[] p = {2, 3}; // f(x) = 3x + 2
        assertEquals(2, Ex1.f(p, 0), Ex1.EPS);
        assertEquals(5, Ex1.f(p, 1), Ex1.EPS);
        assertEquals(-1, Ex1.f(p, -1), Ex1.EPS);
    }

    /**
     * Tests f(x) on a quadratic polynomial.
     */
    @Test
    public void testFQuadratic() {
        double[] p = {1, -2, 1}; // f(x) = x^2 -2x + 1 = (x-1)^2
        assertEquals(0, Ex1.f(p, 1), Ex1.EPS);
        assertEquals(1, Ex1.f(p, 0), Ex1.EPS);
        assertEquals(4, Ex1.f(p, 3), Ex1.EPS);
    }
    /**
     * Tests f(x) on a cubic polynomial.
     */
    @Test
    public void testFCubic() {
        double[] p = {0, 1, 0, -2}; // f(x) = -2x^3 + x
        assertEquals(0, Ex1.f(p, 0), Ex1.EPS);
        assertEquals(-2 + 1, Ex1.f(p, 1), Ex1.EPS);
        assertEquals(16 - 2, Ex1.f(p, -2), Ex1.EPS);
    }

    /**
     * Tests f(x) on the zero polynomial.
     */
    @Test
    public void testFZeroPolynomial() {
        double[] p = {0};
        assertEquals(0, Ex1.f(p, 0), Ex1.EPS);
        assertEquals(0, Ex1.f(p, 15), Ex1.EPS);
        assertEquals(0, Ex1.f(p, -5), Ex1.EPS);
    }

    /**
     * Tests f(x) with a large value of x.
     */
    @Test
    public void testFLargeX() {
        double[] p = {1, 0, 1}; // f(x) = x^2 + 1
        assertEquals(10001, Ex1.f(p, 100), Ex1.EPS);
    }

    /**
     * Tests constructing a polynomial from exactly 2 points.
     */
    @Test
    public void test_polynomFrom2Points() {
        double[] xx = {0, 1};
        double[] yy = {2, 4}; // y = 2x + 2
        double[] p = Ex1.PolynomFromPoints(xx, yy);

        assertEquals(2, p[0], Ex1.EPS);
        assertEquals(2, p[1], Ex1.EPS);
    }

    /**
     * Tests constructing a quadratic polynomial from 3 points forming a parabola.
     */
    @Test
    public void test_polynomFrom3Points_parabola() {
        double[] xx = {0, 1, 2};
        double[] yy = {1, 2, 5}; // y = x^2 + 1
        double[] p = Ex1.PolynomFromPoints(xx, yy);

        assertEquals(1, p[0], Ex1.EPS);
        assertEquals(0, p[1], Ex1.EPS);
        assertEquals(1, p[2], Ex1.EPS);
    }

    /**
     * Tests invalid input for PolynomFromPoints where denominator becomes zero.
     */
    @Test
    public void test_polynomFromPoints_invalid() {
        double[] x = {1, 1};
        double[] y = {2, 2};
        assertNull(Ex1.PolynomFromPoints(x, y));
    }

    /**
     * Tests equals() for two identical polynomials.
     */
    @Test
    public void test_equals_identical() {
        double[] p1 = {1, 2, 3};
        double[] p2 = {1, 2, 3};
        assertTrue(Ex1.equals(p1, p2));
    }

    /**
     * Tests equals() where polynomials have same shape but scaled — should fail.
     */
    @Test
    public void test_equals_scaled_notEqual() {
        double[] p1 = {1, 2, 3};
        double[] p2 = {2, 4, 6};
        assertFalse(Ex1.equals(p1, p2));
    }

    /**
     * Tests equals() where one polynomial has trailing zeros.
     */
    @Test
    public void test_equals_differentSizesButEqual() {
        assertTrue(Ex1.equals(new double[]{1,2,0,0}, new double[]{1,2}));
    }

    /**
     * Tests equals() with values differing by less than EPS.
     */
    @Test
    public void test_equals_withEPS() {
        double[] p1 = {1, 2.0005};
        double[] p2 = {1, 2.0004};
        assertTrue(Ex1.equals(p1, p2));
    }

    /**
     * Tests poly() for a general polynomial with mixed signs.
     */
    @Test
    public void test_poly_basic() {
        double[] p = {2, 0, 3, -1.2};
        String s = Ex1.poly(p);
        assertEquals("-1.2x^3 +3.0x^2 +2.0", s);
    }

    /**
     * Tests poly() on zero polynomial.
     */
    @Test
    public void test_poly_zero() {
        assertEquals("", Ex1.poly(new double[]{0}));
    }

    /**
     * Tests poly() on a simple linear polynomial.
     */
    @Test
    public void test_poly_linear() {
        assertEquals("2.0x +1.0", Ex1.poly(new double[]{1, 2}));
    }

    /**
     * Tests sameValue() for simple linear vs constant — one intersection.
     */
    @Test
    public void test_sameValue_linear() {
        double[] p1 = {1,1};  // x + 1
        double[] p2 = {3};    // 3
        double x = Ex1.sameValue(p1, p2, -10, 10, Ex1.EPS);
        assertEquals(2, x, 0.01);
    }
    /**
     * Tests sameValue() when no intersection exists in the given range.
     */
    @Test
    public void test_sameValue_noSolution() {
        double[] p1 = {1,1};
        double[] p2 = {100};
        assertEquals(-1, Ex1.sameValue(p1, p2, -2, 2, Ex1.EPS));
    }

    /**
     * Tests sameValue() when both functions are constant and never meet.
     */
    @Test
    public void test_sameValue_noCrossing() {
        double[] p1 = {5};
        double[] p2 = {1};
        assertEquals(-1, Ex1.sameValue(p1,p2,-10,10,Ex1.EPS));
    }

    /**
     * Tests sameValue() when x1 > x2 (range is reversed).
     */
    @Test
    public void test_sameValue_reverseRange() {
        double[] p1 = {1,1}; // x+1
        double[] p2 = {3};
        double x = Ex1.sameValue(p1,p2,5,-5,Ex1.EPS);
        assertEquals(2, x, 0.01);
    }

    /**
     * Tests length() for zero-length interval (x1 == x2).
     */
    @Test
    public void test_length_zeroInterval() {
        double[] p = {3,2};
        assertEquals(0, Ex1.length(p, 5, 5, 5));
    }

    /**
     * Tests length() with a single segment approximation.
     */
    @Test
    public void test_length_singleSegment() {
        double[] p = {0,1}; // y=x
        double len = Ex1.length(p, 0, 3, 1);
        assertEquals(Math.sqrt(18), len, Ex1.EPS);
    }

    /**
     * Tests length() on a straight line with many segments.
     */
    @Test
    public void test_length_straightLine() {
        double[] p = {0, 1}; // y = x
        double len = Ex1.length(p, 0, 3, 1000);
        assertEquals(Math.sqrt(2) * 3, len, Ex1.EPS);
    }

    /**
     * Tests length() on a flat constant function.
     */
    @Test
    public void test_length_flat() {
        double[] p = {5};
        assertEquals(3, Ex1.length(p, 0, 3, 10), Ex1.EPS);
    }

    /**
     * Tests area() when x1 == x2 (zero width).
     */
    @Test
    public void test_area_zeroWidth() {
        double[] p1 = {1};
        double[] p2 = {3};
        assertEquals(0, Ex1.area(p1, p2, 5,5,10));
    }
    /**
     * Tests area() when both polynomials are identical — area must be zero.
     */
    @Test
    public void test_area_sameFunction() {
        double[] p = {1,2,3};
        assertEquals(0, Ex1.area(p,p,-10,10,100));
    }

    /**
     * Tests area() for a simple constant gap between two functions.
     */
    @Test
    public void test_area_simple() {
        double[] p1 = {0};     // y = 0
        double[] p2 = {2};     // y = 2
        double a = Ex1.area(p1, p2, 0, 3, 1000);
        assertEquals(6, a, Ex1.EPS);
    }

    /**
     * Tests area() when curves cross inside the interval.
     */
    @Test
    public void test_area_crossing() {
        double[] p1 = {0,1};    // y = x
        double[] p2 = {2};      // y = 2
        double a = Ex1.area(p1, p2, 0, 4, 2000);
        assertTrue(a > 0);
    }

    /**
     * Tests getPolynomFromString() on a full quadratic expression.
     */
    @Test
    public void test_string_basic() {
        String s = "-1.2x^2 +3x +2";
        double[] p = Ex1.getPolynomFromString(s);

        assertEquals(2, p[0], Ex1.EPS);
        assertEquals(3, p[1], Ex1.EPS);
        assertEquals(-1.2, p[2], Ex1.EPS);
    }

    /**
     * Tests getPolynomFromString() with inconsistent spacing.
     */
    @Test
    public void test_string_withSpaces() {
        double[] p = Ex1.getPolynomFromString("  4x  - 2 ");
        assertEquals(-2, p[0], Ex1.EPS);
        assertEquals(4, p[1], Ex1.EPS);
    }

    /**
     * Tests parsing the string "0".
     */
    @Test
    public void test_string_zero() {
        double[] p = Ex1.getPolynomFromString("0");
        assertArrayEquals(new double[]{0}, p);
    }

    /**
     * Tests add() when one or both inputs are null.
     */
    @Test
    public void test_add_null() {
        assertArrayEquals(Ex1.ZERO, Ex1.add(null, null));
        assertArrayEquals(new double[]{1,2}, Ex1.add(new double[]{1,2}, null));
        assertArrayEquals(new double[]{1,2}, Ex1.add(null, new double[]{1,2}));
    }

    /**
     * Tests add() with large positive and negative values cancelling out.
     */
    @Test
    public void test_add_negativeAndLarge() {
        double[] p1 = {1000000, -5};
        double[] p2 = {-1000000, 5};
        assertArrayEquals(new double[]{0,0}, Ex1.add(p1,p2));
    }

    /**
     * Tests add() for two normal polynomials.
     */
    @Test
    public void test_add() {
        double[] p1 = {1,2};
        double[] p2 = {3,4};
        double[] r = Ex1.add(p1,p2);

        assertArrayEquals(new double[]{4,6}, r);
    }
    /**
     * Tests add() where the second polynomial is shorter.
     */
    @Test
    public void test_add_diffSizes() {
        double[] p1 = {1,2,3};
        double[] p2 = {5};
        double[] r = Ex1.add(p1,p2);

        assertArrayEquals(new double[]{6,2,3}, r);
    }

    /**
     * Tests add() with the ZERO polynomial — result must be the same as p1.
     */
    @Test
    public void test_add_zero() {
        double[] p1 = {1,2,3};
        double[] p2 = Ex1.ZERO;

        assertArrayEquals(p1, Ex1.add(p1,p2));
    }

    /**
     * Tests basic subtraction of two equal-length polynomials.
     */
    @Test
    public void test_sub() {
        double[] p1 = {5,4};
        double[] p2 = {2,1};
        double[] r = Ex1.sub(p1,p2);

        assertArrayEquals(new double[]{3,3}, r);
    }

    /**
     * Tests subtraction when the second polynomial is shorter.
     */
    @Test
    public void test_sub_diffSizes() {
        double[] p1 = {10,5,2};
        double[] p2 = {3};
        double[] r = Ex1.sub(p1,p2);

        assertArrayEquals(new double[]{7,5,2}, r);
    }

    /**
     * Tests subtracting ZERO — should return the original polynomial.
     */
    @Test
    public void test_sub_zero() {
        double[] p = {2,3,4};
        assertArrayEquals(p, Ex1.sub(p, Ex1.ZERO));
    }

    /**
     * Tests ZERO - p — should return the negation of p.
     */
    @Test
    public void test_sub_reverseZero() {
        double[] p = {3,5};
        assertArrayEquals(new double[]{-3,-5}, Ex1.sub(Ex1.ZERO, p));
    }

    /**
     * Tests multiplication with NULL inputs — must return null.
     */
    @Test
    public void test_mul_null() {
        assertNull(Ex1.mul(null, new double[]{1}));
        assertNull(Ex1.mul(new double[]{1}, null));
    }

    /**
     * Tests multiplication of a very high-degree polynomial (e.g., x^99).
     */
    @Test
    public void test_mul_largeDegree() {
        double[] p1 = new double[100];
        p1[99] = 1;  // x^99

        double[] p2 = {1,1}; // x + 1

        double[] r = Ex1.mul(p1, p2);

        assertEquals(101, r.length);  // degrees 0..100
        assertEquals(1, r[99]);   // x^99
        assertEquals(1, r[100]);  // x^100
    }

    /**
     * Tests multiplication of (x+1)*(x+1).
     */
    @Test
    public void test_mul() {
        double[] p1 = {1,1}; // x+1
        double[] p2 = {1,1}; // x+1
        double[] r = Ex1.mul(p1,p2); // x^2 +2x +1

        assertArrayEquals(new double[]{1,2,1}, r);
    }

    /**
     * Tests multiplication with ZERO — result must be ZERO.
     */
    @Test
    public void test_mul_zero() {
        double[] r = Ex1.mul(new double[]{3,2}, Ex1.ZERO);
        assertArrayEquals(Ex1.ZERO, r);
    }

    /**
     * Tests general polynomial multiplication with different degrees.
     */
    @Test
    public void test_mul_general() {
        double[] p1 = {2,1};    // x + 2
        double[] p2 = {3,4,5};  // 5x^2 + 4x + 3
        double[] r = Ex1.mul(p1,p2);

        assertArrayEquals(new double[]{6,11,14,5}, r);
    }

    /**
     * Tests derivative of a quadratic polynomial.
     */
    @Test
    public void test_derivative() {
        double[] p = {1,2,3}; // 3x^2 +2x +1
        double[] d = Ex1.derivative(p);
        assertArrayEquals(new double[]{2,6}, d); // f' = 2 + 6x
    }

    /**
     * Tests derivative of a constant polynomial — result must be ZERO.
     */
    @Test
    public void test_derivative_constant() {
        double[] p = {7};
        assertArrayEquals(Ex1.ZERO, Ex1.derivative(p));
    }

    /**
     * Tests normalize(): removing trailing zero coefficients.
     */
    @Test
    public void test_normalize_removeTrailingZeros() {
        double[] p = {1,2,3,0,0};
        assertArrayEquals(new double[]{1,2,3}, Ex1.normalize(p));
    }

    /**
     * Tests normalize() on an all-zero polynomial.
     */
    @Test
    public void test_normalize_allZero() {
        double[] p = {0,0,0};
        assertArrayEquals(Ex1.ZERO, Ex1.normalize(p));
    }

    /**
     * Tests normalize() with null input — should return ZERO.
     */
    @Test
    public void test_normalize_null() {
        assertArrayEquals(Ex1.ZERO, Ex1.normalize(null));
    }

    /**
     * Tests normalize() with many zeros inside.
     */
    @Test
    public void test_normalize_manyZeros() {
        double[] p = {5,0,0,0,0,0};
        assertArrayEquals(new double[]{5}, Ex1.normalize(p));
    }

    /**
     * Tests findRootInSegment() when no root exists in the range.
     */
    @Test
    public void test_findRootInSegment_noRoot() {
        double[] h = {2,1}; // x+2 (positive throughout)
        assertEquals(-1, Ex1.findRootInSegment(h,0,5));
    }

    /**
     * Tests findRootInSegment() when a valid root exists in the interval.
     */
    @Test
    public void test_findRootInSegment() {
        double[] h = {-1,1}; // x-1
        assertEquals(1, Ex1.findRootInSegment(h,0,5), Ex1.EPS);
    }

    /**
     * Tests trapezoid area calculation with negative values.
     */
    @Test
    public void test_trapezoid_negative() {
        assertEquals(5, Ex1.trapezoid(-2, -3, 2), Ex1.EPS); // | -2 | + | -3 |
    }

    /**
     * Tests trapezoid area with positive values.
     */
    @Test
    public void test_trapezoid() {
        assertEquals(5, Ex1.trapezoid(2,3,2), Ex1.EPS);
    }

    /**
     * Tests splitTrapezoid() when the segment crosses zero.
     */
    @Test
    public void test_splitTrapezoid() {
        double area = Ex1.splitTrapezoid(2, -3, 0, 5, 2);
        double expected = 2 + 4.5; // total = 6.5
        assertEquals(expected, area, Ex1.EPS);
    }

    /**
     * Tests splitTrapezoid() when the zero-crossing point is exactly on the edge.
     */
    @Test
    public void test_splitTrapezoid_rootAtEdge() {
        assertEquals(0, Ex1.splitTrapezoid(0,5,2,2,2), Ex1.EPS);
    }
}
