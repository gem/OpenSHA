package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.apache.commons.lang.ArrayUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.imr.attenRelImpl.ToroEtAl_2002_AttenRel;
import org.opensha.sha.imr.attenRelImpl.ToroEtAl_2002_SHARE_AttenRel;
import org.opensha.sha.imr.attenRelImpl.constants.AdjustFactorsSHARE;
import org.opensha.sha.imr.attenRelImpl.constants.ToroEtAl2002Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;

/**
 * Class providing methods for testing {@link ToroEtAl_2002_SHARE_AttenRel}.
 * Tables provided by the original authors. The comparison is done between the
 * tables for {@link ToroEtAl_2002_AttenRel} and the values from
 * {@link ToroEtAl_2002_SHARE_AttenRel} removing the correction.
 */
public class ToroEtAl_2002_SHARE_test implements ParameterChangeWarningListener {

	/** ToroEtAl_2002_AttenRel GMPE (attenuation relationship) */
	private ToroEtAl_2002_SHARE_AttenRel toro2002SHARE = null;

	/**
	 * Table for total standard deviation validation.
	 */
	private static final String SIGMA_TOTAL_HARD_ROCK_TABLE = "Toro02_SIGMAT.txt";

	/**
	 * Table for median ground motion validation. Hard rock median.
	 */
	private static final String MEDIAN_HARD_ROCK_TABLE = "Toro02_MEDIAN.OUT";

	/** Header for meadian tables. */
	private static String[] TABLE_HEADER_MEDIAN = new String[1];

	/** Header for standard deviation tables. */
	private static String[] TABLE_HEADER_STD = new String[1];

	/** Number of columns in test tables for standard deviation. */
	private static final int TABLE_NUM_COL_STD = 12;

	/** Number of columns in test tables for median ground motion value. */
	private static final int TABLE_NUM_COL_MEDIAN = 12;

	/** Number of rows in interface test table. */
	private static final int TABLE_NUM_ROWS = 70;

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
		toro2002SHARE = new ToroEtAl_2002_SHARE_AttenRel(this);
		toro2002SHARE.setParamDefaults();

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
		toro2002SHARE = null;
		medianHardRockTable = null;
		stdTotalTable = null;
	}

	@Test
	public void checkMedianEventOnHardRock() {
		double rake = -90.0;
		validateMedian(rake, medianHardRockTable);
	}

	@Test
	public void checkStdTotal() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_TOTAL, stdTotalTable);
	}

	@Test
	public void checkSetPeriodIndex(){
		double mag = 5.0;
		double rJB = 10.0;
		double rake = 0.0;
		double period = 1.0;
		int iper = ArrayUtils.indexOf(ToroEtAl2002Constants.PERIOD,period);
		toro2002SHARE.getParameter(MagParam.NAME).setValue(mag);
		toro2002SHARE.getParameter(DistanceJBParameter.NAME).setValue(rJB);
		toro2002SHARE.getParameter(RakeParam.NAME).setValue(rake);
		toro2002SHARE.getParameter(PeriodParam.NAME).setValue(period);
		toro2002SHARE.setIntensityMeasure(SA_Param.NAME);
		assertTrue(toro2002SHARE.getMean()==
			toro2002SHARE.getMean(iper, mag, rJB, rake));
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
				double computedMedian = Math.exp(toro2002SHARE.getMean(iper,
						mag, rJB, rake));
				computedMedian = computedMedian
						/ (AdjustFactorsSHARE.AFrock_TORO2002[iper] * toro2002SHARE
								.computeStyleOfFaultingTerm(iper, rake)[2]);
				assertEquals(expectedMedian, computedMedian, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedMedian = table[j][columnDescr.length - 1];
			double computedMedian = Math.exp(toro2002SHARE.getMean(0, mag, rJB,
					rake));
			computedMedian = computedMedian
					/ (AdjustFactorsSHARE.AFrock_TORO2002[0] * toro2002SHARE
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
				double rJB = table[j][1];
				double expectedStd = table[j][i];
				double computedStd = toro2002SHARE.getStdDev(i - 2, mag, rJB,
						stdDevType) / AdjustFactorsSHARE.sig_AFrock_TORO2002[i - 2];
				assertEquals(expectedStd, computedStd, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedStd = table[j][columnDescr.length - 1];
			double computedStd = toro2002SHARE.getStdDev(0, mag, rJB,
					stdDevType) / AdjustFactorsSHARE.sig_AFrock_TORO2002[0];
			assertEquals(expectedStd, computedStd, TOLERANCE);
		}

	}

	@Override
	public void parameterChangeWarning(ParameterChangeWarningEvent event) {
	}
}
