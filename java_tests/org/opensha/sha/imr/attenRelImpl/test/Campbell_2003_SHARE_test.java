package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.*;

import java.io.File;

import org.apache.commons.lang.ArrayUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.param.ParameterAPI;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.attenRelImpl.Campbell_2003_SHARE_AttenRel;
import org.opensha.sha.imr.attenRelImpl.constants.AdjustFactorsSHARE;
import org.opensha.sha.imr.attenRelImpl.constants.Campbell2003Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;

/**
 * Class providing methods for testing {@link Campbell_2003_SHARE_AttenRel}.
 * Tables provided by the original author.
 */
public class Campbell_2003_SHARE_test implements ParameterChangeWarningListener {

	/** Campbell_2003_SHARE_AttenRel GMPE (attenuation relationship) */
	private Campbell_2003_SHARE_AttenRel ca03AttenRel = null;
	/**
	 * Table for total standard deviation validation.
	 */
	private static final String SIGMA_TOTAL_HARD_ROCK_TABLE = "Ca03_SIGMcorr.txt";

	/**
	 * Table for median ground motion validation. Hard rock median.
	 */
	private static final String MEDIAN_HARD_ROCK_TABLE = "Ca03_MEDIAN.OUT";

	/** Header for meadian tables. */
	private static String[] TABLE_HEADER_MEDIAN = new String[1];

	/** Header for standard deviation tables. */
	private static String[] TABLE_HEADER_STD = new String[1];

	/** Number of columns in test tables for standard deviation. */
	private static final int TABLE_NUM_COL_STD = 19;

	/** Number of columns in test tables for median ground motion value. */
	private static final int TABLE_NUM_COL_MEDIAN = 19;

	/** Number of rows in interface test table. */
	private static final int TABLE_NUM_ROWS = 63;

	/** Inter event standard deviation verification table. */
	private static double[][] stdTotalTable = null;

	/** Median ground motion verification table. Normal event on rock. */
	private static double[][] medianHardRockTable = null;

	private static final double TOLERANCE = 1e-3;

	/**
	 * Set up attenuation relationship object, and tables for tests.
	 * 
	 * @throws Exception
	 */
	@Before
	public final void setUp() throws Exception {
		ca03AttenRel = new Campbell_2003_SHARE_AttenRel(this);
		ca03AttenRel.setParamDefaults();

		stdTotalTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_STD];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_TOTAL_HARD_ROCK_TABLE).toURI()),
				stdTotalTable, TABLE_HEADER_STD);
		medianHardRockTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(MEDIAN_HARD_ROCK_TABLE).toURI()),
				medianHardRockTable, TABLE_HEADER_MEDIAN);

	}

	/**
	 * Clean up.
	 */
	@After
	public final void tearDown() {
		ca03AttenRel = null;
		stdTotalTable = null;
		medianHardRockTable = null;
	}

	@Test
	public void checkStdTotal() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_TOTAL, stdTotalTable);
	}

	@Test
	public void checkMedianEventOnHardRock() {
		double rake = -90.0;
		validateMedian(rake, medianHardRockTable);
	}

	@Test
	public void checkSetPeriodIndexGetMean(){
		double mag = 5.0;
		double rRup = 10.0;
		double rake = 0.0;
		double period = 1.0;
		int iper = ArrayUtils.indexOf(Campbell2003Constants.PERIOD, period);
		ca03AttenRel.getParameter(MagParam.NAME).setValue(mag);
		ca03AttenRel.getParameter(DistanceRupParameter.NAME).setValue(rRup);
		ca03AttenRel.getParameter(RakeParam.NAME).setValue(rake);
		ca03AttenRel.getParameter(PeriodParam.NAME).setValue(period);
		ca03AttenRel.setIntensityMeasure(SA_Param.NAME);
		assertTrue(ca03AttenRel.getMean()==
			ca03AttenRel.getMean(iper, mag, rRup, rake));
	}
	
	@Test
	public void checkSetPeriodIndexGetStd(){
		double mag = 5.0;
		double rRup = 10.0;
		double rake = 0.0;
		double period = 1.0;
		String stdDevType = StdDevTypeParam.STD_DEV_TYPE_TOTAL;
		int iper = ArrayUtils.indexOf(Campbell2003Constants.PERIOD, period);
		ca03AttenRel.getParameter(MagParam.NAME).setValue(mag);
		ca03AttenRel.getParameter(DistanceRupParameter.NAME).setValue(rRup);
		ca03AttenRel.getParameter(RakeParam.NAME).setValue(rake);
		ca03AttenRel.getParameter(PeriodParam.NAME).setValue(period);
		ca03AttenRel.setIntensityMeasure(SA_Param.NAME);
		ca03AttenRel.getParameter(StdDevTypeParam.NAME).setValue(stdDevType);
		assertTrue(ca03AttenRel.getStdDev()==
			ca03AttenRel.getStdDev(iper, mag, stdDevType));
	}

	private void validateMedian(double rake, double[][] table) {
		String[] columnDescr = TABLE_HEADER_MEDIAN[0].trim().split("\\s+");
		// check for SA
		for (int i = 3; i < columnDescr.length - 1; i++) {
			for (int j = 0; j < table.length; j++) {
				int iper = i - 2;
				double mag = table[j][0];
				double rJB = table[j][1];
				double expectedMedian = table[j][i];
				double computedMedian = Math.exp(ca03AttenRel.getMean(iper,
						mag, rJB, rake))
						/ (AdjustFactorsSHARE.AFrock_CAMPBELL2003[i-2] * ca03AttenRel
								.computeStyleOfFaultingTerm(iper, rake)[2]);
				assertEquals(expectedMedian, computedMedian, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedMedian = table[j][columnDescr.length - 1];
			double computedMedian = Math.exp(ca03AttenRel.getMean(0, mag, rJB,
					rake))
					/ (AdjustFactorsSHARE.AFrock_CAMPBELL2003[0] * ca03AttenRel
							.computeStyleOfFaultingTerm(0, rake)[2]);
			assertEquals(expectedMedian, computedMedian, TOLERANCE);
		}
	}

	private void validateStdDev(String stdDevType, double[][] table) {
		String[] columnDescr = TABLE_HEADER_STD[0].trim().split("\\s+");
		// check for SA
		for (int i = 3; i < columnDescr.length - 1; i++) {
			for (int j = 0; j < table.length; j++) {
				double mag = table[j][0];
				double expectedStd = table[j][i];
				double computedStd = ca03AttenRel.getStdDev(i - 2, mag,
						stdDevType) / AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003[i - 2];
				assertEquals(expectedStd, computedStd, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double expectedStd = table[j][columnDescr.length - 1];
			double computedStd = ca03AttenRel.getStdDev(0, mag, stdDevType)
					/ AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003[0];
			assertEquals(expectedStd, computedStd, TOLERANCE);
		}
	}

	@Override
	public void parameterChangeWarning(ParameterChangeWarningEvent event) {
	}

    /**
     * This test was added to address
     * https://bugs.launchpad.net/openquake/+bug/942484.
     *
     * The wrong variable was being passed to getStdDev(int, double, String).
     * The 'rake' was being passed instead of 'mag'.
     */
    @Test
    public void testGetStdDev() {

        // Test data was sampled from test_data/Ca03_SIGMcorr.txt
        double mag = 5.0;
        final double [] expectedStdevs = {
                0.6, 0.6, 0.6, 0.623,
                0.633, 0.64, 0.649, 0.658,
                0.662, 0.686, 0.702, 0.7135,
                0.7135, 0.714, 0.7215, 0.731};

        // First, set the stddev type and mag values:
        ParameterAPI magParam = ca03AttenRel.getParameter(MagParam.NAME);
        ParameterChangeEvent magParamChange = new ParameterChangeEvent(
                magParam, MagParam.NAME, null, mag);

        ParameterAPI stddevParam = ca03AttenRel.getParameter(StdDevTypeParam.NAME);
        ParameterChangeEvent stddevParamChange = new ParameterChangeEvent(
                stddevParam, StdDevTypeParam.NAME, null, StdDevTypeParam.STD_DEV_TYPE_TOTAL);
        ca03AttenRel.parameterChange(magParamChange);
        ca03AttenRel.parameterChange(stddevParamChange);

        // Sanity check: make sure
        // expecteStddevs and AFrock_CAMPBELL2003 (periods) have the same length.
        assertEquals(expectedStdevs.length, AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003.length);

       for (int iper = 0; iper < AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003.length; iper++) {
            ca03AttenRel.setIper(iper);
            assertEquals(expectedStdevs[iper], ca03AttenRel.getStdDev() / AdjustFactorsSHARE.sig_AFrock_CAMBPELL2003[iper], 0.0009);
        }
    }
}
