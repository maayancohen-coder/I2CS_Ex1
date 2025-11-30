/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function Computes a polynomial (degree 1 or 2) that passes through 2 or 3 given points.
     * If the input contains 2 points, a line is returned.
     * If it contains 3 points, a quadratic polynomial is computed using a closed formula.
     * If invalid input or parallel x-values, returns null.
     *
     * @param xx array of x-coordinates
     * @param yy array of y-coordinates
     * @return polynomial coefficients in increasing order or null
     *
     * Pseudo code:
     *     if input invalid: return null
     *     if exactly 2 points:
     *         compute slope m = (y2 - y1) / (x2 - x1)
     *         compute intercept c = y1 - m*x1
     *         return [c, m]
     *     if exactly 3 points:
     *         extract x1, y1, x2, y2, x3, y3
     *         compute denominator = (x1-x2)*(x1-x3)*(x2-x3)
     *         if denominator near zero: return null
     *         compute A using formula
     *         compute B using formula
     *         compute C using formula
     *         return [C, B, A]
     */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		    if (lx == 2){
                double x1 = xx[0], y1 = yy[0],
                        x2 = xx[1], y2 = yy[1];
                if (Math.abs(x1 - x2) < EPS) {
                    return null;
                }
                double m = (y2 - y1) / (x2 - x1);
                double c = y1 - m * x1;
                ans = new double[]{c,m};
            }
            else if (lx == 3) {
                double x1 = xx[0], y1 = yy[0],
                        x2 = xx[1], y2 = yy[1],
                        x3 = xx[2], y3 = yy[2];
                double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
                if (Math.abs(denom) < EPS) {
                    return null;
                }
                double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
                double B = (x3 * x3 * (y1 - y2) +
                        x2 * x2 * (y3 - y1) +
                        x1 * x1 * (y2 - y3)) / denom;
                double C = (x2 * x3 * (x2 - x3) * y1 +
                        x3 * x1 * (x3 - x1) * y2 +
                        x1 * x2 * (x1 - x2) * y3) / denom;
                ans = new double[]{C,B,A};
            }
		}
		return ans;
	}
    /**
     * Determines whether two polynomial functions represent the same mathematical function.
     * Two polynomials are considered equal if their f(x) values are equal (up to EPS)
     * for n+1 sample points, where n is the maximum degree of the two polynomials.
     *
     * The method evaluates both polynomials at the points x = 0, 1, 2, ..., n,
     * and checks whether the absolute difference is below EPS for all points.
     *
     * @param p1 the first polynomial, represented as an array of coefficients
     * @param p2 the second polynomial, represented as an array of coefficients
     * @return true if p1 and p2 represent the same polynomial function, false otherwise
     *
     * pseudo code:
     *
     * if both p1 and p2 are null:
     *     return true
     *
     * if exactly one is null:
     *     return false
     *
     * p1 = normalize(p1)
     * p2 = normalize(p2)
     *
     * n = max degree among p1 and p2
     *
     * for x from 0 to n:
     *     compute y1 = f(p1, x)
     *     compute y2 = f(p2, x)
     *     if |y1 - y2| > EPS:
     *         return false
     *
     * return true
     */
	public static boolean equals(double[] p1, double[] p2) {
        if (p1 == null && p2 == null){
            return true;
        }
        if (p1 == null || p2 == null){
            return false;
        }

        p1 = normalize(p1);
        p2 = normalize(p2);

        int n = Math.max(p1.length - 1, p2.length - 1);

        for (int i = 0; i <= n; i++) {
            double v1 = f(p1, i);
            double v2 = f(p2, i);
            if (Math.abs(v1 - v2) > EPS) {
                return false;
            }
        }
        return true;
    }

    /**
     * Converts a polynomial array into a readable String in descending powers.
     * Does not print zero coefficients.
     * Example: input {2,0,3.1,-1.2} → "-1.2x^3 +3.1x^2 +2.0"
     *
     * @param poly the polynomial coefficients
     * @return formatted String of the polynomial
     *
     * Pseudo code:
     *     normalize polynomial
     *     set ans = empty string
     *     loop i from highest degree down to 0:
     *         if coefficient is zero: continue
     *         if ans is empty:
     *             append coefficient
     *         else:
     *             append + or - sign and coefficient
     *         if i > 1: append "x^i"
     *         else if i == 1: append "x"
     *     return ans
     */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            poly = normalize(poly);
            double c;
            for(int i= poly.length-1 ; i>=0 ; i--) {
                c = poly[i];
                if (c == 0.0) {
                    continue;
                }
                if (ans.length() == 0) {
                    ans += c;
                } else {
                    if (c > 0) {
                        ans += " +" + c;
                    } else {
                        ans += " " + c;
                    }
                }
                if (i > 1) {
                    ans += "x^" + i;
                } else if (i == 1) {
                    ans += "x";
                }
            }
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     * Pseudo code:
     *     if either polynomial is null: return -1
     *     compute f1 = p1(x1) - p2(x1)
     *     compute f2 = p1(x2) - p2(x2)
     *     if f1 and f2 have same sign: return -1
     *     while (x2 - x1) > eps:
     *         compute mid = (x1 + x2) / 2
     *         compute fm = p1(mid) - p2(mid)
     *         if |fm| < eps: return mid
     *         if f1 and fm have opposite signs:
     *             x2 = mid, f2 = fm
     *         else:
     *             x1 = mid, f1 = fm
     *     return midpoint of the final interval
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {

        if (p1 == null || p2 == null){
            return -1;
        }

        double f1 = f(p1, x1) - f(p2, x1);
        double f2 = f(p1, x2) - f(p2, x2);

        if (f1 * f2 > 0){
            return -1;
        }

        while (x2 - x1 > eps) {
            double mid = (x1 + x2) / 2;
            double fm = f(p1, mid) - f(p2, mid);

            if (Math.abs(fm) < eps){
                return mid;
            }

            if (f1 * fm <= 0) {
                x2 = mid;
                f2 = fm;
            } else {
                x1 = mid;
                f1 = fm;
            }
        }
        return (x1+x2)/2;
    }

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
     * Pseudo code:
     *     compute dx = x2 - x1
     *     compute delta = dx / numberOfSegments
     *     set ans = 0
     *     for t from x1 to x2 step delta:
     *         compute fx1 = p(t)
     *         compute fx2 = p(t + delta)
     *         compute dy = fx2 - fx1
     *         compute segment_length = sqrt(delta^2 + dy^2)
     *         add to ans
     *     return ans
     */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0;
        double dx = x2 - x1;
        double delta = dx/numberOfSegments;
        double fx1, fx2, dy, dd;
        for (double t = x1; t <= x2; t += delta) {
            fx1 = Ex1.f(p, t);
            fx2 = Ex1.f(p, t+delta);
            dy = fx2 - fx1;
            dd = delta*delta + dy*dy ;
            ans += Math.sqrt(dd);
        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     * Pseudo code:
     *     if invalid input: return 0
     *     if x2 < x1: swap x1 and x2
     *     compute h = p1 - p2 (polynomial subtraction)
     *     compute dx = (x2 - x1) / numberOfTrapezoid
     *     set total = 0
     *     for i from 0 to numberOfTrapezoid - 1:
     *         compute xl = x1 + i*dx
     *         compute xr = xl + dx
     *         compute yL = h(xl)
     *         compute yR = h(xr)
     *         compute root = findRootInSegment(h, xl, xr)
     *         if root == -1:
     *             add trapezoid(|yL|, |yR|, dx) to total
     *         else:
     *             add splitTrapezoid(yL, yR, xl, xr, root) to total
     *     return total
	 */


    public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        if (p1 == null || p2 == null || numberOfTrapezoid <= 0){
            return 0;
        }

        if (x2 < x1) {
            double t = x1;
            x1 = x2;
            x2 = t;
        }

        double[] h = sub(p1, p2);
        double dx = (x2 - x1) / numberOfTrapezoid;
        double total = 0;

        for (int i = 0; i < numberOfTrapezoid; i++) {

            double xl = x1 + i * dx;
            double xr = xl + dx;

            double yL = f(h, xl);
            double yR = f(h, xr);

            double root = findRootInSegment(h, xl, xr);

            if (root == -1) {
                total += trapezoid(yL, yR, dx);
            } else {
                total += splitTrapezoid(yL, yR, xl, xr, root);
            }
        }
        return total;
    }
    /**
     * Detects whether the polynomial h(x) has a root inside [xl, xr].
     * If h(xl) and h(xr) have the same sign, returns -1.
     * Otherwise, returns an approximate root using root_rec.
     *
     * @param h polynomial representing p1 - p2
     * @param xl left x-bound
     * @param xr right x-bound
     * @return root in [xl,xr] or -1 if no sign change
     *
     * Pseudo code:
     *     compute yl = h(xl)
     *     compute yr = h(xr)
     *     if yl*yr > 0:
     *         return -1
     *     return root_rec(h, xl, xr, EPS)
     */
    private static double findRootInSegment(double[] h, double xl, double xr) {
        double yl = f(h, xl);
        double yr = f(h, xr);

        if (yl * yr > 0){
            return -1;
        }

        return root_rec(h, xl, xr, EPS);
    }
    /**
     * Computes the area of a trapezoid with heights y1 and y2 and width dx.
     * Uses |y1| and |y2| because area is always positive.
     *
     * @param y1 value at left endpoint
     * @param y2 value at right endpoint
     * @param dx width of the interval
     * @return trapezoid area
     *
     * Pseudo code:
     *     return (|y1| + |y2|) * dx / 2
     */
    private static double trapezoid(double y1, double y2, double dx) {
        return (Math.abs(y1) + Math.abs(y2)) * dx / 2.0;
    }
    /**
     * Computes area inside a trapezoid where h(x) crosses zero.
     * The interval is split at 'root' into two triangles.
     *
     * @param yL value of h(x) at xl
     * @param yR value of h(x) at xr
     * @param xl left bound
     * @param xr right bound
     * @param root crossing point
     * @return sum of two triangle areas
     *
     * Pseudo code:
     *     compute dx1 = root - xl
     *     compute dx2 = xr - root
     *     compute a1 = |yL|
     *     compute a2 = |yR|
     *     area1 = a1 * dx1 / 2
     *     area2 = a2 * dx2 / 2
     *     return area1 + area2
     */
    private static double splitTrapezoid(double yL, double yR, double xl, double xr, double root) {
        double dx1 = root - xl;
        double dx2 = xr - root;

        double a1 = Math.abs(yL);
        double a2 = Math.abs(yR);

        double area1 = (a1 + 0) * dx1 / 2.0;
        double area2 = (0 + a2) * dx2 / 2.0;

        return area1 + area2;
    }
    /**
     * Converts a polynomial string into an array of coefficients.
     * The returned array represents the polynomial in ascending-power form,
     * meaning index i stores the coefficient of x^i.
     *
     * Example input: "-1.0x^2 +3.0x +2.0"
     * Example output: {2.0, 3.0, -1.0}
     *
     * The method:
     *   - Cleans spaces and ensures every term begins with + or -
     *   - Finds the maximum power to determine array size
     *   - Iterates through the string, extracting:
     *       sign, numeric coefficient, and power
     *   - Accumulates coefficients in the correct index
     *
     * @param p String representing a polynomial
     * @return a double[] representing the polynomial
     *
     * Pseudo code:
     *     clean string using cleanPolyString
     *     if cleaned string invalid: return ZERO
     *     determine max power using getMaxPower
     *     create array of size maxPower+1
     *     set i = 0
     *     while i < string length:
     *         sign = parseSign(string, i)
     *         i++
     *         (coef, newIndex) = parseCoefficient(string, i, sign)
     *         set i = newIndex
     *         power = parsePower(string, i)
     *         i = skipPower(string, i)
     *         ans[power] += coef
     *     normalize array
     *     return array
     */

    public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        String s = cleanPolyString(p);
        if (s == null){
            return ans;
        }

        int maxPower = getMaxPower(s);
        ans = new double[maxPower + 1];

        int i = 0;
        while (i < s.length()) {
            int sign = parseSign(s, i);
            i++;
            double[] coefData = parseCoefficient(s, i, sign);
            double coef = coefData[0];
            i = (int) coefData[1];
            int power = parsePower(s, i);
            i = skipPower(s, i);
            ans[power] += coef;
        }
        normalize(ans);
		return ans;
	}
    /**
     * Removes all spaces from a polynomial string and ensures
     * that the string begins with either '+' or '-'.
     *
     * @param p original string
     * @return cleaned string or null if empty
     *
     * Pseudo code:
     *     if p == null: return null
     *     remove all spaces
     *     if empty: return null
     *     if first char is neither '+' nor '-':
     *         prepend '+'
     *     return cleaned string
     */
    private static String cleanPolyString(String p) {
        if (p == null) return null;
        String s = p.replace(" ", "");
        if (s.length() == 0) return null;

        if (s.charAt(0) != '+' && s.charAt(0) != '-') {
            s = "+" + s;
        }
        return s;
    }
    /**
     * Scans the polynomial string and returns the highest power of x found.
     * If no x is present, returns 0.
     *
     * @param p polynomial string
     * @return maximum exponent
     *
     * Pseudo code:
     *     if p null: return 0
     *     remove spaces
     *     if empty: return 0
     *     maxPower = 0
     *     for each character index i:
     *         if char == 'x':
     *             default power = 1
     *             if next char is '^':
     *                 read digits until non-digit
     *                 convert substring to integer power
     *             update maxPower if necessary
     *     return maxPower
     */
    private static int getMaxPower(String p) {
        if (p == null) {
            return 0;
        }
        String s = p.replace(" ", "");
        if (s.length() == 0) {
            return 0;
        }
        int start, power,maxPower = 0;
        char c;
        for (int i = 0; i < s.length(); i++) {
            c = s.charAt(i);

            if (c == 'x') {
                power = 1;

                if (i + 1 < s.length() && s.charAt(i + 1) == '^') {
                    int j = i + 2;
                    start = j;
                    while (j < s.length() && Character.isDigit(s.charAt(j))) {
                        j++;
                    }
                    String powStr = s.substring(start, j);
                    power = Integer.parseInt(powStr);
                }

                if (power > maxPower) {
                    maxPower = power;
                }
            }
        }

        return maxPower;
    }
    /**
     * Determines the sign of the next polynomial term based on the
     * character at index i. Assumes p[i] is either '+' or '-'.
     *
     * @param s cleaned polynomial string
     * @param i index of sign character
     * @return +1 for '+', -1 for '-'
     *
     * Pseudo code:
     *     if char at i is '-': return -1
     *     else return +1
     */
    private static int parseSign(String s, int i) {
        if (s.charAt(i) == '-'){
            return -1;
        }
        return 1;
    }
    /**
     * Extracts a numeric coefficient starting at position i,
     * reading digits and decimal points until another symbol appears.
     *
     * @param s polynomial string
     * @param i index where the coefficient begins
     * @param sign +1 or -1
     * @return array {coefficientValue, indexAfterCoefficient}
     *
     * Pseudo code:
     *     start = i
     *     while within bounds AND char is digit OR '.':
     *         advance i
     *     extract substring from start to i
     *     convert to double and multiply by sign
     *     return array {value, i}
     */
    private static double[] parseCoefficient(String s, int i, int sign) {
        int start = i;
        while (i < s.length() &&
                (Character.isDigit(s.charAt(i)) || s.charAt(i) == '.')) {
            i++;
        }
        String numStr = s.substring(start, i);
        double coef = Double.parseDouble(numStr) * sign;
        return new double[]{coef, i};
    }
    /**
     * Extracts the exponent of a term.
     * If the next characters form "x^k", the exponent is k.
     * If only "x" appears, exponent is 1.
     * If no x appears, exponent is 0.
     *
     * @param s polynomial string
     * @param i index pointing at 'x' or another character
     * @return exponent of the term
     *
     * Pseudo code:
     *     if out of bounds OR char at i != 'x': return 0
     *     default power = 1
     *     if next char is '^':
     *         advance past '^'
     *         read digits until non-digit
     *         convert digits to integer
     *     return power
     */
    private static int parsePower(String s, int i) {
        if (i >= s.length() || s.charAt(i) != 'x') {
            return 0;
        }
        int power = 1;
        if (i + 1 < s.length() && s.charAt(i + 1) == '^') {
            int j = i + 2;
            int startPow = j;
            while (j < s.length() && Character.isDigit(s.charAt(j))) {
                j++;
            }
            String powStr = s.substring(startPow, j);
            power = Integer.parseInt(powStr);
        }

        return power;
    }
    /**
     * Advances index i past an entire x or x^k expression.
     *
     * Example:
     *   input: "3x^12", i = indexOf('x')
     *   returns index after the digits '12'
     *
     * @param s polynomial string
     * @param i index pointing at 'x'
     * @return new index after skipping power section
     *
     * Pseudo code:
     *     if out of bounds OR char != 'x': return i
     *     i++  // skip 'x'
     *     if next char == '^':
     *         i++ // skip '^'
     *         while char is digit:
     *             i++
     *     return i
     */
    private static int skipPower(String s, int i) {
        if (i >= s.length() || s.charAt(i) != 'x') {
            return i;
        }
        i++;
        if (i < s.length() && s.charAt(i) == '^') {
            i++;
            while (i < s.length() && Character.isDigit(s.charAt(i))) {
                i++;
            }
        }
        return i;
    }
    /**
     * Computes p1 + p2 and returns the resulting polynomial.
     * The result has length equal to the longer polynomial.
     *
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return polynomial representing p1 + p2
     *
     * Pseudo code:
     *     if both null: return ZERO
     *     determine maxLength and minLength
     *     if p1 shorter than p2: swap them
     *     copy longer polynomial into ans
     *     for i from 0 to minLength-1:
     *         ans[i] += p2[i]
     *     return ans
     */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        if(p1 == null && p2 == null) {
            ans = ZERO;
        }
        int maxLength = Math.max(p1.length,p2.length);
        int minLength = Math.min(p1.length,p2.length);
        double [] temp = new double[maxLength];
        if (p1.length < p2.length) {
            for (int i = 0; i < p2.length; i++) {
                temp[i] = p2[i];
            }
            p2 = p1;
            p1 = temp;
        }
        if (p1.length == maxLength) {
            ans = new double[p1.length];
            for(int i=0;i<p1.length;i++) {
                ans[i] = p1[i];
            }
        }
        else {
            ans = new double[p2.length];
            for(int i=0;i<p2.length;i++) {
                ans[i] = p2[i];
            }
        }
        for (int i =0; i<minLength; i++) {
            ans[i] += p2[i];
        }

		return ans;
	}
    /**
     * Computes the product of two polynomials p1 * p2.
     *
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return resulting polynomial
     *
     * Pseudo code:
     *     if either null: return null
     *     if either is ZERO: return ZERO
     *     create array of size (len1 + len2 - 1)
     *     for each i in p1:
     *         for each j in p2:
     *             ans[i+j] += p1[i] * p2[j]
     *     normalize result
     *     return ans
     */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        if (p1 == null || p2 == null){
            return null;
        }
        if (equals(p1, ZERO) || equals(p2, ZERO)){
            return ZERO;
        }
        int n1 = p1.length;
        int n2 = p2.length;
        ans = new double[n1 + n2 - 1];
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                ans[i + j] += p1[i] * p2[j];
            }
        }
        normalize(ans);
        return ans;
	}
    /**
     * Computes the derivative of polynomial p.
     *
     * Example:
     *   p = {c0, c1, c2, c3}
     *   derivative = {c1*1, c2*2, c3*3}
     *
     * @param po input polynomial
     * @return derivative polynomial
     *
     * Pseudo code:
     *     if null or length == 0: return ZERO
     *     create ans of size (length-1)
     *     for i from 1 to end:
     *         ans[i-1] = po[i] * i
     *     return ans
     */
	public static double[] derivative (double[] po) {
		double [] ans;//
        if (po.length==0 || po == null) {
            ans=ZERO;
        } else {
            ans = new double[po.length-1];
        }
        for (int i = 1 ; i < po.length ; i++) {
            ans [i-1] = po[i]*i;
        }
		return ans;
	}
    /**
     * Removes trailing zero coefficients from a polynomial.
     * Ensures minimal length representation.
     *
     * @param p polynomial
     * @return trimmed polynomial with no trailing zeros
     *
     * Pseudo code:
     *     if null or empty: return ZERO
     *     find last non-zero index
     *     if none found: return ZERO
     *     copy coefficients up to last index
     *     return new array
     */
    private static double[] normalize(double[] p) {
        if (p == null || p.length == 0){
            return ZERO;
        }
        int last = p.length-1;
        while (last > 0 && p[last] == 0.0) {
            last--;
        }

        if (last == 0 && p[last] == 0.0) {
            return ZERO;
        }

        double[] ans = new double[last + 1];
        for (int i = 0; i <= last; i++) ans[i] = p[i];
        return ans;
    }
    /**
     * Computes the polynomial p1 - p2.
     * Implemented by multiplying p2 by -1 and then adding.
     *
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return p1 - p2
     *
     * Pseudo code:
     *     if both null: return ZERO
     *     if p1 null: treat as ZERO
     *     if p2 null: treat as ZERO
     *     negP2 = mul(p2, {-1})
     *     ans = add(p1, negP2)
     *     normalize ans
     *     return ans
     */
    public static double [] sub(double [] p1, double[] p2) {
        if (p1 == null && p2 == null) {
            return ZERO;
        }
        if (p1 == null) {
            p1 = ZERO;
        }
        if (p2 == null) {
            p2 = ZERO;
        }
        // -1 כפול p2
        double[] minusOne = {-1};
        double[] negP2 = mul(p2, minusOne);   // (-1) * p2

        double[] ans = add(p1, negP2);

        return normalize(ans);
    }
}
