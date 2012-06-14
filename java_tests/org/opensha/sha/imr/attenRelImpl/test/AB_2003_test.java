package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.imr.attenRelImpl.AB_2003_AttenRel;
import org.opensha.sha.imr.param.EqkRuptureParams.FocalDepthParam;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.TectonicRegionTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Class providing methods for testing {@link AB_2003_AttenRel}. Tables (for the
 * global model, without corrections for Japan/Cascadia) were provided by Celine
 * Beauval (<celine.beauval@obs.ujf-grenoble.fr>) using matSHA.
 */
public class AB_2003_test implements ParameterChangeWarningListener {

	/** Atkinson and Boore 2003 attenuation relationship. */
	private AB_2003_AttenRel ab2003AttenRel = null;

	/**
	 * Table for peak ground acceleration interface test. 1st column: rupture
	 * distance. 2nd column: PGA (g), for M=8.5, site class B 3rd column: PGA
	 * (g), for M=7.5, site class C 4th column: PGA (g), for M=6.5, site class D
	 * Hypocentral depth = 20 km
	 */
	private static final String PGA_INTERFACE_TABLE_FILE = "AtkinsonBoore2003Global-Interface-PGA-g.dat";

	/**
	 * Table for spectral acceleration interface test. 1st column: rupture
	 * distance. 2nd column: SA-1Hz (g), for M=8.5, site class B 3rd column:
	 * SA-1Hz (g), for M=7.5, site class C 4th column: SA-1Hz (g), for M=6.5,
	 * site class D Hypocentral depth = 20 km
	 */
	private static final String SA_INTERFACE_TABLE_FILE = "AtkinsonBoore2003Global-Interface-1Hz-g.dat";

	/**
	 * Table for peak ground acceleration intraslab test. 1st column: rupture
	 * distance. 2nd column: PGA (g), for M=8.0, site class B 3rd column: PGA
	 * (g), for M=7.0, site class C 4th column: PGA (g), for M=6.0, site class D
	 * Hypocentral depth = 60 km
	 */
	private static final String PGA_INTRASLAB_TABLE_FILE = "AtkinsonBoore2003Global-M876-PGA-INTRASLAB-g.dat";

	/**
	 * Table for spectral acceleration intraslab test. 2nd column: SA-1Hz (g),
	 * for M=8.0, site class B 3rd column: SA-1Hz (g), for M=7.0, site class C
	 * 4th column: SA-1Hz (g), for M=6.0, site class D Hypocentral depth = 60 km
	 */
	private static final String SA_INTRASLAB_TABLE_FILE = "AtkinsonBoore2003Global-M876-1Hz-INTRASLAB-g.dat";

	/** Number of columns in test tables. */
	private static final int TABLE_NUM_COL = 4;

	/** Number of rows in interface test table. */
	private static final int INTERFACE_TABLE_NUM_ROWS = 24;

	/** Number of rows in intraslab test table. */
	private static final int INTRA_SLAB_TABLE_NUM_ROWS = 20;

	/**
	 * Vs30 values for interface/intraSlab tests. Values corresponding to class
	 * B,C, and D respectively.
	 */
	private static final double[] VS_30 = new double[] { 800.0, 500.0, 200.0 };
	/**
	 * Magnitude values for interface tests.
	 */
	private static final double[] INTERFACE_MAGNITUDE_VALUES = new double[] {
			8.5, 7.5, 6.5 };
	/**
	 * Magnitude values for intraslab tests.
	 */
	private static final double[] INTRASLAB_MAGNITUDE_VALUES = new double[] {
			8.0, 7.0, 6.0 };

	/** Peak ground acceleration table for interface test. */
	private static double[][] pgaInterfaceTable = null;

	/** Spectral acceleration table for interface test. */
	private static double[][] saInterfaceTable = null;

	/** Peak ground acceleration table for intraSlab test. */
	private static double[][] pgaIntraSlabTable = null;

	/** Spectral acceleration table for intraSlab test. */
	private static double[][] saIntraSlabTable = null;

	/** Hypocentral depth for interface test. */
	private static final double INTERFACE_HYPO_DEPTH = 20.0;

	/** Hypocentral depth for intraSlab test. */
	private static final double INTRA_SLAB_HYPO_DEPTH = 60.0;

	/** Tolerance level (expressed in percentage). */
	private static final double TOLERANCE = 1e-3;

	/**
	 * Set up attenuation relationship object, and tables for tests.
	 * 
	 * @throws Exception
	 */
	@Before
	public final void setUp() throws Exception {
		ab2003AttenRel = new AB_2003_AttenRel(this);
		ab2003AttenRel.setParamDefaults();
		pgaInterfaceTable = new double[INTERFACE_TABLE_NUM_ROWS][TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(new File(ClassLoader
				.getSystemResource(PGA_INTERFACE_TABLE_FILE).toURI()),
				pgaInterfaceTable);
		saInterfaceTable = new double[INTERFACE_TABLE_NUM_ROWS][TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(new File(ClassLoader
				.getSystemResource(SA_INTERFACE_TABLE_FILE).toURI()),
				saInterfaceTable);
		pgaIntraSlabTable = new double[INTRA_SLAB_TABLE_NUM_ROWS][TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(new File(ClassLoader
				.getSystemResource(PGA_INTRASLAB_TABLE_FILE).toURI()),
				pgaIntraSlabTable);
		saIntraSlabTable = new double[INTRA_SLAB_TABLE_NUM_ROWS][TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(new File(ClassLoader
				.getSystemResource(SA_INTRASLAB_TABLE_FILE).toURI()),
				saIntraSlabTable);
	}

	/**
	 * Clean up.
	 */
	@After
	public final void tearDown() {
		ab2003AttenRel = null;
		pgaInterfaceTable = null;
		saInterfaceTable = null;
		pgaIntraSlabTable = null;
		saIntraSlabTable = null;
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw=8.0, intraslab event,
	 * site type NEHRP B, hypocentral depth 60 km.
	 */
	@Test
	public final void sa1HzMw8IntraSlabNerphBHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 0;
		int siteTypeIndex = 0;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 1;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saIntraSlabTable);
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw7.0, intraslab event, site
	 * type Nerph C, hypocentral depth 60 km.
	 */
	@Test
	public final void sa1HzMw7IntraSlabNerphCHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 1;
		int siteTypeIndex = 1;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 2;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saIntraSlabTable);
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw = 6.0, intraslab event,
	 * site type NEHRP D, hypocentral depth 60 km.
	 */
	@Test
	public final void sa1HzMw65IntraSlabNerphDHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 2;
		int siteTypeIndex = 2;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 3;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saIntraSlabTable);
	}

	/**
	 * Check median pga, for Mw = 8.0, intraslab event, site type NEHRP B,
	 * hypocentral depth 60 km.
	 */
	@Test
	public final void pgaMw8IntraSlabNerphBHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 0;
		int siteTypeIndex = 0;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 1;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaIntraSlabTable);
	}

	/**
	 * Check median pga for Mw = 7.0, intra slab event, site type Nerph C,
	 * hypocentral depth 60 km.
	 */
	@Test
	public final void pgaMw8IntraSlabNerphCHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 1;
		int siteTypeIndex = 1;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 2;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaIntraSlabTable);
	}

	/**
	 * Check median pga for Mw = 6.0, intraslab event, site type NEHRP D,
	 * hypocentral depth 60 km.
	 */
	@Test
	public final void pgaMw6IntraSlabNerphDHypodepth60() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_SLAB
				.toString();
		double hypocentralDepth = INTRA_SLAB_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 2;
		int siteTypeIndex = 2;
		double magnitude = INTRASLAB_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 3;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaIntraSlabTable);
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw=8.5, interface event,
	 * site type NEHRP B, hypocentral depth 20 km.
	 */
	@Test
	public final void sa1HzMw85InterfaceNerphBHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 0;
		int siteTypeIndex = 0;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 1;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saInterfaceTable);
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw7.5, interface event, site
	 * type Nerph C, hypocentral depth 20 km.
	 */
	@Test
	public final void sa1HzMw75InterfaceNerphCHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 1;
		int siteTypeIndex = 1;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 2;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saInterfaceTable);
	}

	/**
	 * Check median spectral acceleration (1Hz) for Mw = 6.5, interface event,
	 * site type NEHRP D, hypocentral depth 20 km.
	 */
	@Test
	public final void sa1HzMw65InterfaceNerphDHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 5;
		int magnitudeIndex = 2;
		int siteTypeIndex = 2;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 3;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, saInterfaceTable);
	}

	/**
	 * Check median pga for Mw=8.5, interface event, site type NEHRP B,
	 * hypocentral depth 20 km.
	 */
	@Test
	public final void pgaMw85InterfaceNerphBHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 0;
		int siteTypeIndex = 0;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 1;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaInterfaceTable);
	}

	/**
	 * Check median pga for Mw7.5, interface event, site type Nerph C,
	 * hypocentral depth 20 km.
	 */
	@Test
	public final void pgaMw75InterfaceNerphCHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 1;
		int siteTypeIndex = 1;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 2;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaInterfaceTable);
	}

	/**
	 * Check median pga for Mw = 6.5, interface event, site type NEHRP D,
	 * hypocentral depth 20 km.
	 */
	@Test
	public final void pgaMw65InterfaceNerphDHypodepth20() {

		String tectonicRegionType = TectonicRegionType.SUBDUCTION_INTERFACE
				.toString();
		double hypocentralDepth = INTERFACE_HYPO_DEPTH;
		int periodIndex = 0;
		int magnitudeIndex = 2;
		int siteTypeIndex = 2;
		double magnitude = INTERFACE_MAGNITUDE_VALUES[magnitudeIndex];
		double vs30 = VS_30[siteTypeIndex];
		int expectedResultIndex = 3;

		validateAgainstTable(tectonicRegionType, hypocentralDepth, periodIndex,
				magnitude, vs30, expectedResultIndex, pgaInterfaceTable);
	}
	
	/**
	 * Check SA at T = 4 s
	 */
	@Test
	public final void sa4s(){
		ab2003AttenRel.setIntensityMeasure(SA_Param.NAME);
		ab2003AttenRel.getParameter(MagParam.NAME).setValue(6.0);
		ab2003AttenRel.getParameter(DistanceRupParameter.NAME).setValue(10.0);
		ab2003AttenRel.getParameter(Vs30_Param.NAME).setValue(800.0);
		ab2003AttenRel.getParameter(TectonicRegionTypeParam.NAME).
			setValue(TectonicRegionType.SUBDUCTION_INTERFACE.toString());
		ab2003AttenRel.getParameter(FocalDepthParam.NAME).setValue(20.0);
		
		ab2003AttenRel.getParameter(PeriodParam.NAME).setValue(3.0);
		double mean3s = ab2003AttenRel.getMean();
		
		ab2003AttenRel.getParameter(PeriodParam.NAME).setValue(4.0);
		double mean4s = ab2003AttenRel.getMean();
		
		assertTrue(mean4s  == mean3s / 0.550);
	}

	/**
	 * Compare median ground motion againt values in table.
	 * 
	 * @param tectonicRegionType
	 * @param hypocentralDepth
	 * @param periodIndex
	 * @param magnitudeIndex
	 * @param siteTypeIndex
	 * @param expectedResultIndex
	 * @param table
	 */
	private void validateAgainstTable(final String tectonicRegionType,
			final double hypocentralDepth, final int periodIndex,
			final double mag, final double vs30, final int expectedResultIndex,
			final double[][] table) {
		for (int i = 0; i < table.length; i++) {

			double distance = table[i][0];

			double predicted = Math.exp(ab2003AttenRel.getMean(periodIndex,
					mag, distance, vs30, tectonicRegionType, hypocentralDepth));
			double expected = table[i][expectedResultIndex];
			double percentageDifference = Math.abs((expected - predicted)
					/ expected) * 100;

			String msg = "distance: " + distance + ", magnitude: " + mag
					+ ", vs30: " + vs30 + ", tectonic region type: "
					+ tectonicRegionType + ", hypocentral depth: "
					+ hypocentralDepth + ", expected: " + expected
					+ ", predicted: " + predicted + ",percentage difference: "
					+ percentageDifference;
			assertTrue(msg, percentageDifference < TOLERANCE);
		}
	}

	@Override
	public void parameterChangeWarning(final ParameterChangeWarningEvent event) {
	}

}
