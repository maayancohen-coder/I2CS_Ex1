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
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
     * pseudo code:
     * for
     * denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
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
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
        if (p1 == null && p2 == null){
            return true;
        }
        if (p1 == null || p2 == null){
            return false;
        }

        double temp = 0;
        boolean result = true;
        for (double i = 0; i < p1.length; i++) {
            temp = f(p1,i) - f(p2,i);
            if (Math.abs(temp) > EPS) {
                result = false;

            }
        }
        return result;
        /**if (p1.length != p2.length) return false;

        for (int i = 0; i < p1.length; i++) {
            if (Math.abs(p1[i] - p2[i]) > EPS) {
                return false;
            }
        }
        return true;*/
    }

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
     * לתקן - הפונקציות צריכות פשוט לתת את אותם ערכים בקירוב של אפסילון כדי להיות שוות
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
     * fix!!!!
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans = -1;
        double[] t = {-1};
        double[] h = add(p2, mul(p1, t));
        for (double xi = x1; xi <= x2; xi++) {
            if (f(h, xi) <=EPS) {
                return ans = x1;
            }
        }
        return ans;
    }
		/**double ans = x1;
        double left = x1;
        double right = x2;
        double yLeft = f(p1,left) - f(p2,left);
        double xMid, yMid;
        while (Math.abs(left - right) > eps) {
            xMid = (left + right) / 2.0;
            yMid = f(p1,xMid) - f(p2,xMid);
            if (Math.abs(yMid) < eps) {
                ans = xMid;
                break;
            }
            if (yLeft * yMid <= 0){
                right = xMid;
            } else {
                left = xMid;
                yLeft = yMid;
            }
            ans = (left + right) / 2.0;
        }
		return ans;
	}*/
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
     *
     * double h = (x2 - x1)/numberOfTrapezoid
     * for (double xi=x1 ; xi<=x2 ; xi+= h){
     *     gx1  = Math.abs(f(p2,xi) - f(p1,xi))
     *     gx2 = Math.abs(f(p2,xi+h) - f(p1,xi+h))
     *     ans += ((gx1 + gx2)*h)/2.0
     * }
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        if (numberOfTrapezoid<=0){
            return ans;
        }
        if (x2<x1){
            double temp = x1;
            x1 = x2;
            x2 = temp;
        }
        double h = (x2 - x1)/numberOfTrapezoid;
        double gx1,gx2;

        for (double xi=x1 ; xi<=x2 ; xi+= h){
            gx1  = Math.abs(f(p2,xi) - f(p1,xi));
            gx2 = Math.abs(f(p1,xi+h) - f(p2,xi+h));
            ans += ((gx1 + gx2)*h)/2.0;
        }
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
     * for
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
    private static String cleanPolyString(String p) {
        if (p == null) return null;
        String s = p.replace(" ", "");
        if (s.length() == 0) return null;

        if (s.charAt(0) != '+' && s.charAt(0) != '-') {
            s = "+" + s;
        }
        return s;
    }
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
    private static int parseSign(String s, int i) {
        if (s.charAt(i) == '-'){
            return -1;
        }
        return 1;
    }
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
         * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
         * @param p1
         * @param p2
         * @return
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
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
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
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
     * for (i=0, i<p0.length, i++)
     *  if i = 0
     *  i++
     *  ans [i-1] = po[i]*i
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

    private static double[] normalize(double[] p) {
        if (p == null || p.length == 0){
            return ZERO;
        }
        int last = p.length-1;
        while (last > 0 && Math.abs(p[last]) < EPS) {
            last--;
        }
        if (last == 0 && Math.abs(p[last]) < EPS) {
            return ZERO;
        }
        if (last == p.length) {
            return p;
        }
        double [] ans = new double [last + 1];
        for (int i = 0 ; i <= last ; i++) {
            ans [i] = p[i] ;
        }
        return ans;
    }
}
