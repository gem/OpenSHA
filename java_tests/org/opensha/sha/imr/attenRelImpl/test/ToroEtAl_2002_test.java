package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.imr.attenRelImpl.ToroEtAl_2002_AttenRel;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;

/**
 * Class providing methods for testing {@link ToroEtAl_2002_AttenRel}. Tables
 * provided by the original authors.
 */
public class ToroEtAl_2002_test implements ParameterChangeWarningListener {

	/** ToroEtAl_2002_AttenRel GMPE (attenuation relationship) */
	private ToroEtAl_2002_AttenRel toro2002AtenRel = null;

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
		toro2002AtenRel = new ToroEtAl_2002_AttenRel(this);
		toro2002AtenRel.setParamDefaults();

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
		toro2002AtenRel = null;
		medianHardRockTable = null;
		stdTotalTable = null;
	}

	@Test
	public void checkMedianEventOnHardRock() {
		validateMedian(medianHardRockTable);
	}

	@Test
	public void checkStdTotal() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_TOTAL, stdTotalTable);
	}
	
	/**
	 * Check SA at T = 3 and 4 s
	 */
	@Test
	public final void sa3_4s(){
		toro2002AtenRel.setIntensityMeasure(SA_Param.NAME);
		toro2002AtenRel.getParameter(MagParam.NAME).setValue(6.0);
		toro2002AtenRel.getParameter(DistanceJBParameter.NAME).setValue(10.0);
		
		toro2002AtenRel.getParameter(PeriodParam.NAME).setValue(2.0);
		double mean2s = toro2002AtenRel.getMean();
		
		toro2002AtenRel.getParameter(PeriodParam.NAME).setValue(3.0);
		double mean3s = toro2002AtenRel.getMean();
		
		toro2002AtenRel.getParameter(PeriodParam.NAME).setValue(4.0);
		double mean4s = toro2002AtenRel.getMean();
		
		assertTrue(mean3s  == mean2s / 0.385);
		assertTrue(mean4s  == mean2s / 0.1483);
	}

	private void validateMedian(double[][] table) {
		String[] columnDescr = TABLE_HEADER_MEDIAN[0].trim().split("\\s+");
		// check for SA
		for (int i = 3; i < columnDescr.length - 1; i++) {
			for (int j = 0; j < table.length; j++) {
				int iper = i-2;
				double mag = table[j][0];
				double rJB = table[j][1];
				double expectedMedian = table[j][i];
				double computedMedian = Math.exp(toro2002AtenRel.getMean(iper, mag, rJB));
				assertEquals(expectedMedian, computedMedian, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedMedian = table[j][columnDescr.length - 1];
			double computedMedian = Math.exp(toro2002AtenRel.getMean(0, mag, rJB));
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
				double computedStd = toro2002AtenRel.getStdDev(i-2, mag, 
						rJB, stdDevType);
				assertEquals(expectedStd, computedStd, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedStd = table[j][columnDescr.length - 1];
			double computedStd = toro2002AtenRel.getStdDev(0, mag, 
					rJB, stdDevType);
			assertEquals(expectedStd, computedStd, TOLERANCE);
		}
	}

	@Override
	public void parameterChangeWarning(ParameterChangeWarningEvent event) {
	}
}
