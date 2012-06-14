package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.imr.attenRelImpl.AkB_2010_AttenRel;
import org.opensha.sha.imr.attenRelImpl.constants.AkB2010Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.EqkRuptureParams.RakeParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceJBParameter;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;

/**
 * Class providing methods for testing {@link AkB_2010_AttenRel}. Tables
 * provided by the original authors.
 */
public class AkB_2010_test implements ParameterChangeWarningListener {

	/** Akkar and Bommer 2010 attenuation relationship. */
	private AkB_2010_AttenRel akb2010AttenRel = null;

	/**
	 * Table for inter-event standard deviation validation.
	 */
	private static final String SIGMA_INTER_NM_ROCK_TABLE = "AB10_SIGMAINTER_NM_ROCK.OUT";

	/**
	 * Table for intra-event standard deviation validation.
	 */
	private static final String SIGMA_INTRA_NM_ROCK_TABLE = "AB10_SIGMAINTRA_NM_ROCK.OUT";

	/**
	 * Table for total standard deviation validation.
	 */
	private static final String SIGMA_TOTAL_NM_ROCK_TABLE = "AB10_SIGMAT_NM_ROCK.OUT";

	/**
	 * Table for median ground motion validation. Normal event on rock.
	 */
	private static final String MEDIAN_NM_ROCK_TABLE = "AB2010_MEDIAN_NM_ROCK.OUT";
	
	/**
	 * Table for median ground motion validation. Normal event on soft soil.
	 */
	private static final String SIGMA_MEDIAN_NM_SOFTSOIL_TABLE = "AB2010_MEDIAN_NM_SOFTSOIL.OUT";
	
	/**
	 * Table for median ground motion validation. Normal event on stiff soil.
	 */
	private static final String SIGMA_MEDIAN_NM_STIFFSOIL_TABLE = "AB2010_MEDIAN_NM_STIFFSOIL.OUT";
	
	/**
	 * Table for median ground motion validation. Reverse event on rock.
	 */
	private static final String SIGMA_MEDIAN_REVERSE_ROCK_TABLE = "AB2010_MEDIAN_RR_ROCK.OUT";
	
	/**
	 * Table for median ground motion validation. Reverse event on soft soil.
	 */
	private static final String SIGMA_MEDIAN_REVERSE_SOFTSOIL_TABLE = "AB2010_MEDIAN_RR_SOFTSOIL.OUT";
	
	/**
	 * Table for median ground motion validation. Reverse event on stiff soil.
	 */
	private static final String SIGMA_MEDIAN_REVERSE_STIFFSOIL_TABLE = "AB2010_MEDIAN_RR_STIFFSOIL.OUT";
	
	/**
	 * Table for median ground motion validation. Strike-slip event on rock.
	 */
	private static final String SIGMA_MEDIAN_STRIKESLIP_ROCK_TABLE = "AB2010_MEDIAN_SS_ROCK.OUT";
	
	/**
	 * Table for median ground motion validation. Strike-slip event on soft soil.
	 */
	private static final String SIGMA_MEDIAN_STRIKESLIP_SOFTSOIL_TABLE = "AB2010_MEDIAN_SS_SOFTSOIL.OUT";
	
	/**
	 * Table for median ground motion validation. Strike-slip event on stiff soil.
	 */
	private static final String SIGMA_MEDIAN_STRIKESLIP_STIFFSOIL_TABLE = "AB2010_MEDIAN_SS_STIFFSOIL.OUT";

	/** Header for meadian tables. */
	private static String[] TABLE_HEADER_MEDIAN = new String[1];
	
	/** Header for standard deviation tables. */
	private static String[] TABLE_HEADER_STD = new String[1];

	/** Number of columns in test tables for standard deviation. */
	private static final int TABLE_NUM_COL_STD = 69;

	/** Number of columns in test tables for median ground motion value. */
	private static final int TABLE_NUM_COL_MEDIAN = 68;

	/** Number of rows in interface test table. */
	private static final int TABLE_NUM_ROWS = 42;

	/** Inter event standard deviation verification table. */
	private static double[][] stdInterTable = null;

	/** Inter event standard deviation verification table. */
	private static double[][] stdIntraTable = null;

	/** Inter event standard deviation verification table. */
	private static double[][] stdTotalTable = null;

	/** Median ground motion verification table. Normal event on rock. */
	private static double[][] medianNormalRockTable = null;
	
	/** Median ground motion verification table. Normal event on soft soil. */
	private static double[][] medianNormalSoftSoilTable = null;
	
	/** Median ground motion verification table. Normal event on stiff soil. */
	private static double[][] medianNormalStiffSoilTable = null;
	
	/** Median ground motion verification table. Reverse event on rock. */
	private static double[][] medianReverseRockTable = null;
	
	/** Median ground motion verification table. Reverse event on soft soil. */
	private static double[][] medianReverseSoftSoilTable = null;
	
	/** Median ground motion verification table. Reverse event on stiff soil. */
	private static double[][] medianReverseStiffSoilTable = null;
	
	/** Median ground motion verification table. Strike slip event on rock. */
	private static double[][] medianStrikeSlipRockTable = null;
	
	/** Median ground motion verification table. Strike slip event on soft soil. */
	private static double[][] medianStrikeSlipSoftSoilTable = null;
	
	/** Median ground motion verification table. Strike slip event on stiff soil. */
	private static double[][] medianStrikeSlipStiffSoilTable = null;

	private static final double TOLERANCE = 1e-4;

	/**
	 * Set up attenuation relationship object, and tables for tests.
	 * 
	 * @throws Exception
	 */
	@Before
	public final void setUp() throws Exception {
		akb2010AttenRel = new AkB_2010_AttenRel(this);
		akb2010AttenRel.setParamDefaults();
		stdInterTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_STD];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_INTER_NM_ROCK_TABLE).toURI()),
				stdInterTable, TABLE_HEADER_STD);
		stdIntraTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_STD];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_INTRA_NM_ROCK_TABLE).toURI()),
				stdIntraTable, TABLE_HEADER_STD);
		stdTotalTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_STD];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_TOTAL_NM_ROCK_TABLE).toURI()),
				stdTotalTable, TABLE_HEADER_STD);
		medianNormalRockTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(MEDIAN_NM_ROCK_TABLE).toURI()),
				medianNormalRockTable, TABLE_HEADER_MEDIAN);
		medianNormalSoftSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_NM_SOFTSOIL_TABLE).toURI()),
				medianNormalSoftSoilTable, TABLE_HEADER_MEDIAN);
		medianNormalStiffSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_NM_STIFFSOIL_TABLE).toURI()),
				medianNormalStiffSoilTable, TABLE_HEADER_MEDIAN);
		medianReverseRockTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_REVERSE_ROCK_TABLE).toURI()),
				medianReverseRockTable, TABLE_HEADER_MEDIAN);
		medianReverseSoftSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_REVERSE_SOFTSOIL_TABLE).toURI()),
				medianReverseSoftSoilTable, TABLE_HEADER_MEDIAN);
		medianReverseStiffSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_REVERSE_STIFFSOIL_TABLE).toURI()),
				medianReverseStiffSoilTable, TABLE_HEADER_MEDIAN);
		medianStrikeSlipRockTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_STRIKESLIP_ROCK_TABLE).toURI()),
				medianStrikeSlipRockTable, TABLE_HEADER_MEDIAN);
		medianStrikeSlipSoftSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_STRIKESLIP_SOFTSOIL_TABLE).toURI()),
				medianStrikeSlipSoftSoilTable, TABLE_HEADER_MEDIAN);
		medianStrikeSlipStiffSoilTable = new double[TABLE_NUM_ROWS][TABLE_NUM_COL_MEDIAN];
		AttenRelTestHelper.readNumericTableWithHeader(new File(ClassLoader
				.getSystemResource(SIGMA_MEDIAN_STRIKESLIP_STIFFSOIL_TABLE).toURI()),
				medianStrikeSlipStiffSoilTable, TABLE_HEADER_MEDIAN);
	}

	/**
	 * Clean up.
	 */
	@After
	public final void tearDown() {
		akb2010AttenRel = null;
		stdInterTable = null;
		stdIntraTable = null;
		stdTotalTable = null;
		medianNormalRockTable = null;
		medianNormalSoftSoilTable = null;
		medianNormalStiffSoilTable = null;
		medianReverseRockTable = null;
		medianReverseSoftSoilTable = null;
		medianStrikeSlipRockTable = null;
		medianStrikeSlipSoftSoilTable = null;
		medianStrikeSlipStiffSoilTable = null;
	}

	@Test
	public void checkStdInter() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_INTER, stdInterTable);
	}

	@Test
	public void checkStdIntra() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_INTRA, stdIntraTable);
	}

	@Test
	public void checkStdTotal() {
		validateStdDev(StdDevTypeParam.STD_DEV_TYPE_TOTAL, stdTotalTable);
	}

	@Test
	public void checkMedianNormalEventOnRock() {
		double vs30 = 800.0;
		double rake = -90.0;
		validateMedian(vs30, rake, medianNormalRockTable);
	}
	
	@Test
	public void checkMedianNormalEventOnSoftSoil() {
		double vs30 = 100.0;
		double rake = -90.0;
		validateMedian(vs30, rake, medianNormalSoftSoilTable);
	}
	
	@Test
	public void checkMedianNormalEventOnStiffSoil() {
		double vs30 = 400.0;
		double rake = -90.0;
		validateMedian(vs30, rake, medianNormalStiffSoilTable);
	}
	
	@Test
	public void checkMedianReverseEventOnRock() {
		double vs30 = 800.0;
		double rake = 90.0;
		validateMedian(vs30, rake, medianReverseRockTable);
	}
	
	@Test
	public void checkMedianReverseEventOnSoftSoil() {
		double vs30 = 100.0;
		double rake = 90.0;
		validateMedian(vs30, rake, medianReverseSoftSoilTable);
	}
	
	@Test
	public void checkMedianReverseEventOnStiffSoil() {
		double vs30 = 400.0;
		double rake = 90.0;
		validateMedian(vs30, rake, medianReverseStiffSoilTable);
	}
	
	@Test
	public void checkMedianStrikeSlipEventOnRock() {
		double vs30 = 800.0;
		double rake = 0.0;
		validateMedian(vs30, rake, medianStrikeSlipRockTable);
	}
	
	@Test
	public void checkMedianStrikeSlipEventOnSoftSoil() {
		double vs30 = 100.0;
		double rake = 0.0;
		validateMedian(vs30, rake, medianStrikeSlipSoftSoilTable);
	}
	
	@Test
	public void checkMedianStrikeSlipEventOnStiffSoil() {
		double vs30 = 400.0;
		double rake = 0.0;
		validateMedian(vs30, rake, medianStrikeSlipStiffSoilTable);
	}
	
	/**
	 * Check SA at T = 4 s
	 */
	@Test
	public final void sa4s(){
		akb2010AttenRel.setIntensityMeasure(SA_Param.NAME);
		akb2010AttenRel.getParameter(MagParam.NAME).setValue(6.0);
		akb2010AttenRel.getParameter(RakeParam.NAME).setValue(0.0);
		akb2010AttenRel.getParameter(DistanceJBParameter.NAME).setValue(10.0);
		akb2010AttenRel.getParameter(Vs30_Param.NAME).setValue(800.0);
		
		akb2010AttenRel.getParameter(PeriodParam.NAME).setValue(3.0);
		double mean3s = akb2010AttenRel.getMean();
		
		akb2010AttenRel.getParameter(PeriodParam.NAME).setValue(4.0);
		double mean4s = akb2010AttenRel.getMean();
		
		assertTrue(mean4s  == mean3s / 0.80);
	}

	private void validateMedian(double vs30, double rake, double[][] table) {
		String[] columnDescr = TABLE_HEADER_MEDIAN[0].split(" ");
		// check for SA
		for (int i = 2; i < columnDescr.length - 2; i++) {
			for (int j = 0; j < table.length; j++) {
				int iper = i;
				double mag = table[j][0];
				double rJB = table[j][1];
				double expectedMedian = table[j][i];
				double computedMedian = Math.exp(akb2010AttenRel.getMean(iper,
						mag, rJB, vs30, rake));
				assertEquals(expectedMedian, computedMedian, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedMedian = table[j][columnDescr.length - 2];
			double computedMedian = Math.exp(akb2010AttenRel.getMean(1, mag,
					rJB, vs30, rake));
			assertEquals(expectedMedian, computedMedian, TOLERANCE);
		}
		// check for PGV
		for (int j = 0; j < table.length; j++) {
			double mag = table[j][0];
			double rJB = table[j][1];
			double expectedMedian = table[j][columnDescr.length - 1];
			double computedMedian = Math.exp(akb2010AttenRel.getMean(0, mag,
					rJB, vs30, rake));
			assertEquals(expectedMedian, computedMedian, TOLERANCE);
		}
	}

	private void validateStdDev(String stdDevType, double[][] table) {
		String[] columnDescr = TABLE_HEADER_STD[0].split(" ");
		// check for SA
		for (int i = 3; i < columnDescr.length - 2; i++) {
			for (int j = 0; j < table.length; j++) {
				double expectedStd = table[j][i];
				double computedStd = akb2010AttenRel.getStdDev(i - 1,
						stdDevType) / AkB2010Constants.LOG10_2_LN;
				assertEquals(expectedStd, computedStd, TOLERANCE);
			}
		}
		// check for PGA
		for (int j = 0; j < table.length; j++) {
			double expectedStd = table[j][columnDescr.length - 2];
			double computedStd = akb2010AttenRel.getStdDev(1, stdDevType)
					/ AkB2010Constants.LOG10_2_LN;
			assertEquals(expectedStd, computedStd, TOLERANCE);
		}
		// check for PGV
		for (int j = 0; j < table.length; j++) {
			double expectedStd = table[j][columnDescr.length - 1];
			double computedStd = akb2010AttenRel.getStdDev(0, stdDevType)
					/ AkB2010Constants.LOG10_2_LN;
			assertEquals(expectedStd, computedStd, TOLERANCE);
		}
	}

	@Override
	public void parameterChangeWarning(ParameterChangeWarningEvent event) {
	}
}
