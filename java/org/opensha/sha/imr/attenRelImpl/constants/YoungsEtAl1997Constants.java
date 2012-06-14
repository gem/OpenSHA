package org.opensha.sha.imr.attenRelImpl.constants;

public class YoungsEtAl1997Constants {

	public static final double[] PERIOD_SOIL = {  
		    0.00, 0.075, 0.10, 0.20, 0.30, 0.40, 0.500, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00 };
	public static final double[] C1_SOIL = { -0.438, 
		    0.000, 2.400, 2.516, 1.549, 0.793, 0.144, -0.438, -1.704, -2.870, -5.101, 
		    -6.433, -6.672, -7.618 };
	public static final double[] C2_SOIL = {  
		    -0.0019, -0.0019, -0.0019, -0.0020, -0.0020, -0.0035, -0.0048, -0.0066, -0.0114, -0.0164,
			-0.0221, -0.0235, -0.0235 };
	public static final double[] C3_SOIL = {  
		    -2.329, -2.697, -2.697, -2.464, -2.327, -2.230, -2.140, -1.952, -1.785, -1.470, 
		    -1.290, -1.347, -1.272};
	public static final double[] C4_SOIL = {  
		    1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.50,
		    1.55, 1.65, 1.65 };
	public static final double[] C5_SOIL = { 
		    -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10,
		    -0.10, -0.10, -0.10 };

	public static final double[] PERIOD_ROCK = { 
		    0.00, 0.075, 0.10, 0.20, 0.30, 0.40, 0.500, 0.75, 1.00, 1.50, 
		    2.00, 3.00, 4.00 };
	
	public static final double[] C1_ROCK = {  
		     0.000, 1.275, 1.188, 0.722, 0.246, -0.115, -0.400, -1.149, -1.736, -2.634,
		     -3.328, -4.511, -4.511 };
	
	public static final double[] C2_ROCK = {  
		    0.0000, 0.0000, -0.0011, -0.0027, -0.0036, -0.0043, -0.0048, -0.0057, -0.0064, -0.0073, 
		   -0.0080,-0.0089,-0.0089};
	
	public static final double[] C3_ROCK = { 
		    -2.552, -2.707, -2.655, -2.528, -2.454, -2.401, -2.360, -2.286, -2.234, -2.160, 
		    -2.107, -2.033, -2.033 };
	
	public static final double[] C4_ROCK = { 
		    1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.45, 1.50, 
		    1.55, 1.65, 1.65};
	
	public static final double[] C5_ROCK = { 
		   -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, -0.10, 
		   -0.10, -0.10, -0.10 };

	public static final double A1_ROCK = 0.2418;
	public static final double A2_ROCK = 1.4140;
	public static final double A3_ROCK = 10.000;
	public static final double A4_ROCK = 1.7818;
	public static final double A5_ROCK = 0.554;
	public static final double A6_ROCK = 0.00607;
	public static final double A7_ROCK = 0.3846;

	public static final double A1_SOIL = -0.6687;
	public static final double A2_SOIL = 1.438;
	public static final double A3_SOIL = 10.000;
	public static final double A4_SOIL = 1.097;
	public static final double A5_SOIL = 0.617;
	public static final double A6_SOIL = 0.00648;
	public static final double A7_SOIL = 0.3643;

	public static final Double MAG_WARN_MIN = new Double(5.00);
	public static final Double MAG_WARN_MAX = new Double(8.20);

	public static final Double DISTANCE_RUP_WARN_MIN = new Double(8.5);
	public static final Double DISTANCE_RUP_WARN_MAX = new Double(551.0);
	/**
	 * Class (B+C) Vs30 >= 360 
	 * Class (D+E) Vs30 < 360 
	 */
	public static final double SOIL_TYPE_SOFT_UPPER_BOUND = 360.0; 

	/**
	 * cm/s/s to g conversion factor.
	 */
	public static final double CMS_TO_G_CONVERSION_FACTOR = 1.0/981.0;
	
    /**
	* Decay factor to extend the spectral periods to 4sec (for rock)
	*/
	public static final double T3sec_TO_T4sec_DECAYFACTOR = 0.399;
}
