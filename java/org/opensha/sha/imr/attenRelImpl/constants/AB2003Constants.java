package org.opensha.sha.imr.attenRelImpl.constants;

public class AB2003Constants {

	/**
	 * Supported frequency values (Hz).
	 */
	public static final double[] FREQ = {Double.POSITIVE_INFINITY,
		25.0000, 10.0000, 5.00000, 2.50000, 1.00000, 0.50000, 0.33000, 0.25000};
	/**
	 * Supported period values (s).
	 */
	public static final double[] PERIOD = {0.00000, 0.04000, 0.10000,
		0.20000, 0.40000, 1.00000, 2.00000, 3.00000, 4.00000};
	/**
	 * Interface c1 coefficients.
	 */
	public static final double[] INTER_C1 = {2.99100, 2.87530, 2.77890,
		2.66380, 2.52490, 2.14420, 2.19070, 2.30100, 2.30100};
	/**
	 * Interface c2 coefficients.
	 */
	public static final double[] INTER_C2 = {0.03525, 0.07052, 0.09841,
		0.12386, 0.14770, 0.13450, 0.07148, 0.02237, 0.02237};
	/**
	 * Interface c3 coefficients.
	 */
	public static final double[] INTER_C3 = {0.00759, 0.01004, 0.00974,
		0.00884, 0.00728, 0.00521, 0.00224, 0.00012, 0.00012};
	/**
	 * Interface c4 coefficients.
	 */
	public static final double[] INTER_C4 = {-0.00206, -0.00278, -0.00287,
		-0.00280, -0.00235, -0.00110, 0.00000, 0.00000, 0.00000};
	/**
	 * Interface c5 coefficients.
	 */
	public static final double[] INTER_C5 = {0.19000, 0.15000, 0.15000,
		0.15000, 0.13000, 0.10000, 0.10000, 0.10000, 0.10000};
	/**
	 * Interface c6 coefficients.
	 */
	public static final double[] INTER_C6 = {0.24000, 0.20000, 0.23000,
		0.27000, 0.37000, 0.30000, 0.25000, 0.25000, 0.25000};
	/**
	 * Interface c7 coefficients.
	 */
	public static final double[] INTER_C7 = {0.29000, 0.20000, 0.20000,
		0.25000, 0.38000, 0.55000, 0.40000, 0.36000, 0.36000};
	/**
	 * Interface total standard deviation.
	 */
	public static final double[] INTER_TOTAL_STD = {0.23000, 0.26000,
		0.27000, 0.28000, 0.29000, 0.34000, 0.34000, 0.36000, 0.36000};
	/**
	 * Interface intra-event standard deviation.
	 */
	public static final double[] INTER_INTRAEVENT_STD = {0.20000, 0.22000,
		0.25000, 0.25000, 0.25000, 0.28000, 0.29000, 0.31000, 0.31000};
	/**
	 * Interface inter-event standard deviation.
	 */
	public static final double[] INTER_INTEREVENT_STD = {0.11000, 0.14000,
		0.10000, 0.13000, 0.15000, 0.19000, 0.18000, 0.18000, 0.18000};
	/**
	 * Intraslab c1 coefficients.
	 */
	public static final double[] INTRA_C1 = {-0.04713, 0.50697, 0.43928,
		0.51589, 0.00545, -1.02133, -2.39234, -3.70012, -3.70012};
	/**
	 * Intraslab c2 coefficients.
	 */
	public static final double[] INTRA_C2 = {0.69090, 0.63273, 0.66675,
		0.69186, 0.77270, 0.87890, 0.99640, 1.11690, 1.11690};
	/**
	 * Intraslab c3 coefficients.
	 */
	public static final double[] INTRA_C3 = {0.01130, 0.01275, 0.01080,
		0.00572, 0.00173, 0.00130, 0.00364, 0.00615, 0.00615};
	/**
	 * Intraslab c4 coefficients.
	 */
	public static final double[] INTRA_C4 = {-0.00202, -0.00234,
		-0.00219, -0.00192, -0.00178, -0.00173, -0.00118, -0.00045, -0.00045};
	/**
	 * Intraslab c5 coefficients.
	 */
	public static final double[] INTRA_C5 = {0.19000, 0.15000, 0.15000,
		0.15000, 0.13000, 0.10000, 0.10000, 0.10000, 0.10000};
	/**
	 * Intraslab c6 coefficients.
	 */
	public static final double[] INTRA_C6 = {0.24000, 0.20000, 0.23000,
		0.27000, 0.37000, 0.30000, 0.25000, 0.25000, 0.25000};
	/**
	 * Intraslab c7 coefficients.
	 */
	public static final double[] INTRA_C7 = {0.29000, 0.20000, 0.20000,
		0.25000, 0.38000, 0.55000, 0.40000, 0.36000, 0.36000};
	/**
	 * Intraslab total standard deviation.
	 */
	public static final double[] INTRA_TOTAL_STD = {0.27000, 0.25000,
		0.28000, 0.28000, 0.28000, 0.29000, 0.30000, 0.30000, 0.30000};
	/**
	 * Intraslab intra event standard deviation.
	 */
	public static final double[] INTRA_INTRAEVENT_STD = {0.23000, 0.24000,
		0.27000, 0.26000, 0.26000, 0.27000, 0.28000, 0.29000, 0.29000};
	/**
	 * Intraslab inter event standard deviation.
	 */
	public static final double[] INTRA_INTEREVENT_STD = {0.14000, 0.07000,
		0.07000, 0.10000, 0.10000, 0.11000, 0.11000, 0.08000, 0.08000};
	/**
	 * log10 to natural log conversion factor.
	 */
	public static final double LOG10_2_LN = Math.log(10.0);
	/**
	 * Minimum magnitude.
	 */
	public static final Double MAG_WARN_MIN = new Double(5.5);
	/**
	 * Maximum magnitude.
	 */
	public static final Double MAG_WARN_MAX = new Double(8.3);

	/**
	 * Minimum rupture distance.
	 */
	public static final Double DISTANCE_RUP_WARN_MIN = new Double(11.0);
	/**
	 * Maximum rupture distance.
	 */
	public static final Double DISTANCE_RUP_WARN_MAX = new Double(550.0);
	/**
	 * Minimum hypocentral depth.
	 */
	public static final Double DEPTH_HYPO_WARN_MIN = new Double(0.0);
	/**
	 * Maximum hypocentral depth.
	 */
	public static final Double DEPTH_HYPO_WARN_MAX = new Double(125.0);
	
	/**
	 * Threshold hypocentral depth.
	 */
	public static final double THRESHOLD_HYPO_DEPTH = 100.0;
	/**
	 * Threshold magnitude for interafce events.
	 */
	public static final double THRESHOLD_MAG_INTERFACE = 8.5;
	/**
	 * Threshold magnitude for intraslab events.
	 */
	public static final double THRESHOLD_MAG_INTRASLAB = 8.0;
	/**
	 * NEHRP E Vs30 upper bound
	 */
	public static final double NEHRP_E_UPPER_BOUND = 180.0; 
	/**
	 * NEHRP D Vs30 upper bound
	 */
	public static final double NEHRP_D_UPPER_BOUND = 360.0; 
	/**
	 * NEHRP C Vs30 upper bound
	 */
	public static final double NEHRP_C_UPPER_BOUND = 760.0; 
	/**
	 * Interface geometric spreading factor f1
	 */
	public static final double INTERFACE_GEOM_SPREAD_FACTOR1 = 1.2;
	/**
	 * Interface geometric spreading factor f2
	 */
	public static final double INTERFACE_GEOM_SPREAD_FACTOR2 = 0.18;
	/**
	 * Intraslab geometric spreading factor f1
	 */
	public static final double INTRASLAB_GEOM_SPREAD_FACTOR1 = 0.301;
	/**
	 * Intraslab geometric spreading factor f2
	 */
	public static final double INTRASLAB_GEOM_SPREAD_FACTOR2 = 0.01;
	/**
	 * Near source saturation term factor 1
	 */
	public static final double NEAR_SOURCE_SATURATION_FACTOR1 = 0.00724;
	/**
	 * Near source saturation term factor 2
	 */
	public static final double NEAR_SOURCE_SATURATION_FACTOR2 = 0.507;
	/**
	 * cm/s to g conversion factor.
	 */
	public static final double CMS_TO_G_CONVERSION_FACTOR = 1.0/981.0;
	/**
	 * weights for correction
	 */
	public static final double[] CORRECTION_WEIGHTS =
		new double[]{ 0.333, 0.667};
	/**
	 * decay factor to extend spectral period from 3 to 4 s
	 */
	public static final double T3sec_TO_T4sec_DECAYFACTOR = 0.550;
}
