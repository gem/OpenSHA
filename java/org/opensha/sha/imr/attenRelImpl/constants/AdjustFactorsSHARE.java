package org.opensha.sha.imr.attenRelImpl.constants;

public class AdjustFactorsSHARE {

	/**
	 * Rock adjustment coefficients coefficients to correct the hard-rock
	 * soil(vs30>=2800m/sec) to rock (vs30>=800m/sec) obtained from
	 * interpolation of the coefficients presented in table 9 Drouet et al
	 * (2010)
	 */
	public static final double[] AFrock_TORO2002 = { 0.735106, 0.419632,
			0.477379, 0.888509, 1.197291, 1.308267, 1.30118, 1.265762,
			1.215779, 1.215779, 1.215779};

	public static final double[] AFrock_CAMPBELL2003 = { 0.735106, 0.474275,
			0.423049, 0.550323, 0.730061, 0.888509, 1.094622, 1.197291,
			1.288309, 1.311421, 1.298212, 1.265762, 1.197583, 1.215779,
			1.215779, 1.215779 };

	/**
	 * Adjustment factors of the total standard deviation when hard rock (vs30
	 * =2800m/s) to refference rock (vs30 =800m/s). Interpolated from Table 9
	 * (Drouet et al 2010) Note: the adjustment factors for periods of 3 and 4
	 * sec are identical with those for 2sec - SHARE experts decision
	 */
	public static final double[] sig_AFrock_TORO2002 = { 0.338916, 0.289785,
			0.320650, 0.352442, 0.281552, 0.198424, 0.1910745, 0.154327,
			0.155520, 0.155520, 0.155520};

	public static final double[] sig_AFrock_CAMBPELL2003 = { 0.338916,
			0.283461, 0.289785, 0.345375, 0.365490, 0.352442, 0.315477,
			0.281552, 0.231117, 0.175678, 0.149567, 0.154327, 0.191149,
			0.155520, 0.155520, 0.155520 };
	/**
	 * Style-of-faulting adjustment coefficients, obtained from Table 7 (Drouet
	 * et al 2010) and using a cubic spline interpolation
	 */
	public static final double[] Frss_TORO2002 = {1.220000, 0.935198, 0.907936,
			1.080000, 1.190000, 1.230000, 1.230000, 1.177500, 1.140000, 1.140000,
			1.140000};
	
	public static final double[] Frss_CAMPBELL2003 = { 1.220000, 1.080745, 0.986646,
			0.931209, 0.910332, 0.985998, 1.080000, 1.150000, 1.190000,
			1.230000, 1.191702, 1.177500, 1.194705, 1.140000, 1.140000,
			1.140000 };

	public static final double Fnss = 0.95;
	public static final double pN = 0.01;
	public static final double pR = 0.81;

	/**
	 * NORMAL - STYLE of FAULTING (Definition based on rake angle minimum value)
	 */
	public static final double FLT_TYPE_NORMAL_RAKE_LOWER = -120.00;
	/**
	 * NORMAL - STYLE of FAULTING (Definition based on rake angle maximum value)
	 */
	public static final double FLT_TYPE_NORMAL_RAKE_UPPER = -60.00;
	/**
	 * Reverse - STYLE of FAULTING (Definition based on rake angle minimum
	 * value)
	 */
	public static final double FLT_TYPE_REVERSE_RAKE_LOWER = 30.00;
	/**
	 * Reverse - STYLE of FAULTING (Definition based on rake angle maximum
	 * value)
	 */
	public static final double FLT_TYPE_REVERSE_RAKE_UPPER = 150.00;

	public static void main(String[] args) {
		System.out.println(AFrock_TORO2002.length);
		System.out.println(sig_AFrock_TORO2002.length);
		System.out.println(Frss_CAMPBELL2003.length);
	}

}
