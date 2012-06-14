package org.opensha.sha.imr.attenRelImpl.test;

import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.opensha.commons.data.function.ArbitrarilyDiscretizedFunc;
import org.opensha.commons.param.event.ParameterChangeWarningEvent;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.imr.attenRelImpl.YoungsEtAl_1997_AttenRel;
import org.opensha.sha.imr.param.EqkRuptureParams.FocalDepthParam;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.TectonicRegionTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.TectonicRegionType;

/**
 * Class providing methods for testing {@link YoungsEtAl_1997_AttenRel}.
 */
public class YoungsEtAl_1997_test implements ParameterChangeWarningListener {
	
	/** Youngs et al. 1997 GMPE */
	private YoungsEtAl_1997_AttenRel youngsEtAl1997AttenRel= null;
	
	/** Table for peak ground acceleration interface test. 
	 * 1st column: distance (km)
	 * 2nd column: PGA(g), for Mw=8.8
	 * 3rd column: PGA(g), for Mw=8.0
	 * 4th column: PGA(g), for Mw=7.0
	 * tectonic region type is INTERFACE, site class "A",
	 * hypocentral depth 25 km.
	 */
	private static final String PGA_INTERFACE_TABLE_FILE =
		"Youngs-interface-PGA-g.dat";
	
	/** Table for spectral acceleration interface test.
	 * 1st column: distance (km)
	 * 2nd column: SA 1Hz (g), for Mw=8.8
	 * 3rd column: SA 1Hz (g), for Mw=8.0
	 * 4th column: SA 1Hz (g), for Mw=7.0
	 * tectonic region type is INTERFACE, site class "A",
	 * hypocentral depth 25 km. 
	 */
	private static final String SA_INTERFACE_TABLE_FILE =
		"Youngs-interface-1Hz-g.dat";
	
	/** Table for peak ground acceleration intraslab test. 
	 * 1st column: distance (km)
	 * 2nd column: PGA(g), for Mw=8.0
	 * 3rd column: PGA(g), for Mw=7.5
	 * 4th column: PGA(g), for Mw=7.0
	 * tectonic region type is INTRASLAB, site class "A",
	 * hypocentral depth 60 km.
	 */
	private static final String PGA_INTRASLAB_TABLE_FILE =
		"Youngs-intraslab-PGA-g.dat";
	
	/** Table for spectral acceleration intraslab test. 
	 * 1st column: distance (km)
	 * 2nd column: SA 1Hz (g), for Mw=8.0
	 * 3rd column: SA 1Hz (g), for Mw=7.5
	 * 4th column: SA 1Hz (g), for Mw=7.0
	 * tectonic region type is INTRASLAB, site class "A",
	 * hypocentral depth 60 km.
	 */
	private static final String SA_INTRASLAB_TABLE_FILE =
		"Youngs-intraslab-1Hz-g.dat";
	
	/** Number of columns in test tables. */
	private static final int TABLE_NUM_COL = 4;
	
	/** Number of rows in interface test table. */
	private static final int TABLE_NUM_ROWS = 24;
	
	/** Magnitude values for interface tests. */
	private static final double[] INTERFACE_MAGNITUDE_VALUES =
		new double[] { 8.8, 8.0, 7.0 };
	
	/** Magnitude values for intraslab tests. */
	private static final double[] INTRASLAB_MAGNITUDE_VALUES =
		new double[] { 8.0, 7.5, 7.0 };
	
	/** Vs30 corresponding to rock **/
	private static final double VS_30 = 800.0;
	
	/** Peak ground acceleration table for interface test. */
	private static double[][] pgaInterfaceTable = null;
	
	/** Spectral acceleration table for interface test. */
	private static double[][] saInterfaceTable = null;
	
	/** Peak ground acceleration table for intraSlab test. */
	private static double[][] pgaIntraSlabTable = null;
	
	/** Spectral acceleration table for intraSlab test. */
	private static double[][] saIntraSlabTable = null;

	/** Hypocentral depth for interface test. */
	private static final double INTERFACE_HYPO_DEPTH = 25.0;
	
	/** Hypocentral depth for intraSlab test. */
	private static final double INTRA_SLAB_HYPO_DEPTH = 60.0;
	
	/** Tolerance level (expressed in percentage). */
	private static final double TOLERANCE = 0.5;
	
	/**
	 * Set up attenuation relationship object, and tables for tests.
	 * 
	 * @throws Exception
	 */
	@Before
	public final void setUp() throws Exception {
		youngsEtAl1997AttenRel = new YoungsEtAl_1997_AttenRel(this);
		youngsEtAl1997AttenRel.setParamDefaults();
		pgaInterfaceTable = 
			new double[YoungsEtAl_1997_test.TABLE_NUM_ROWS]
			           [YoungsEtAl_1997_test.TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(
				new File(ClassLoader.getSystemResource(
						YoungsEtAl_1997_test.PGA_INTERFACE_TABLE_FILE).toURI()),
				pgaInterfaceTable);
		saInterfaceTable = 
			new double[YoungsEtAl_1997_test.TABLE_NUM_ROWS]
			           [YoungsEtAl_1997_test.TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(
				new File(ClassLoader.getSystemResource(
						YoungsEtAl_1997_test.SA_INTERFACE_TABLE_FILE).toURI()),
				saInterfaceTable);
		pgaIntraSlabTable =
			new double[YoungsEtAl_1997_test.TABLE_NUM_ROWS]
			           [YoungsEtAl_1997_test.TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(
				new File(ClassLoader.getSystemResource(
						YoungsEtAl_1997_test.PGA_INTRASLAB_TABLE_FILE).toURI()),
				pgaIntraSlabTable);
		saIntraSlabTable =
			new double[YoungsEtAl_1997_test.TABLE_NUM_ROWS]
			           [YoungsEtAl_1997_test.TABLE_NUM_COL];
		AttenRelTestHelper.readNumericTable(
				new File(ClassLoader.getSystemResource(
						YoungsEtAl_1997_test.SA_INTRASLAB_TABLE_FILE).toURI()),
				saIntraSlabTable);
	}
	
	/**
	 * Clean up.
	 */
	@After
	public final void tearDown() {
		youngsEtAl1997AttenRel = null;
		pgaInterfaceTable = null;
		saInterfaceTable = null;
		pgaIntraSlabTable = null;
		saIntraSlabTable = null;
	}
	
	/**
	 * Check for PGA from Mw=8.8 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void pgaInterfaceM88HypoDepth25(){
		int iper = 0;
		double mag = INTERFACE_MAGNITUDE_VALUES[0];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 1;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaInterfaceTable);
	}
	
	/**
	 * Check for PGA from Mw=8.0 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void pgaInterfaceM8HypoDepth25(){
		int iper = 0;
		double mag = INTERFACE_MAGNITUDE_VALUES[1];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 2;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaInterfaceTable);
	}
	
	/**
	 * Check for PGA from Mw=7.0 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void pgaInterfaceM7HypoDepth25(){
		int iper = 0;
		double mag = INTERFACE_MAGNITUDE_VALUES[2];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 3;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaInterfaceTable);
	}
	
	/**
	 * Check for SA (1Hz) from Mw=8.8 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void saInterfaceM88HypoDepth25(){
		int iper = 8;
		double mag = INTERFACE_MAGNITUDE_VALUES[0];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 1;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saInterfaceTable);
	}
	
	/**
	 * Check for SA from Mw=8.0 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void saInterfaceM8HypoDepth25(){
		int iper = 8;
		double mag = INTERFACE_MAGNITUDE_VALUES[1];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 2;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saInterfaceTable);
	}
	
	/**
	 * Check for SA from Mw=7.0 for interface events,
	 * hypocentral depth = 25 km.
	 */
	@Test
	public void saInterfaceM7HypoDepth25(){
		int iper = 8;
		double mag = INTERFACE_MAGNITUDE_VALUES[2];
		String tecRegType = TectonicRegionType.SUBDUCTION_INTERFACE.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTERFACE_HYPO_DEPTH;
		
		int expectedResultIndex = 3;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saInterfaceTable);
	}
	
	/**
	 * Check for PGA from Mw=8.0 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void pgaIntraslabM8HypoDepth60(){
		int iper = 0;
		double mag = INTRASLAB_MAGNITUDE_VALUES[0];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 1;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaIntraSlabTable);
	}
	
	/**
	 * Check for PGA from Mw=7.5 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void pgaIntraslabM75HypoDepth60(){
		int iper = 0;
		double mag = INTRASLAB_MAGNITUDE_VALUES[1];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 2;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaIntraSlabTable);
	}
	
	/**
	 * Check for PGA from Mw=7.0 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void pgaIntraslabM7HypoDepth60(){
		int iper = 0;
		double mag = INTRASLAB_MAGNITUDE_VALUES[2];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 3;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, pgaIntraSlabTable);
	}
	
	/**
	 * Check for SA (1Hz) from Mw=8.0 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void saIntraslabM8HypoDepth60(){
		int iper = 8;
		double mag = INTRASLAB_MAGNITUDE_VALUES[0];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 1;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saIntraSlabTable);
	}
	
	/**
	 * Check for SA from Mw=7.5 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void saIntraslabM75HypoDepth60(){
		int iper = 8;
		double mag = INTRASLAB_MAGNITUDE_VALUES[1];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 2;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saIntraSlabTable);
	}
	
	/**
	 * Check for SA from Mw=7.0 for intraslab events,
	 * hypocentral depth = 60 km.
	 */
	@Test
	public void saIntraslabM7HypoDepth60(){
		int iper = 8;
		double mag = INTRASLAB_MAGNITUDE_VALUES[2];
		String tecRegType = TectonicRegionType.SUBDUCTION_SLAB.toString();
		String siteTypeParam = "Rock";
		double hypoDep = INTRA_SLAB_HYPO_DEPTH;
		
		int expectedResultIndex = 3;
		
		validateAgainstTable(tecRegType,
				hypoDep, iper,
				mag, siteTypeParam,
				expectedResultIndex, saIntraSlabTable);
	}
	
	/**
	 * Check SA at T = 4 s
	 */
	@Test
	public final void sa4s(){
		youngsEtAl1997AttenRel.setIntensityMeasure(SA_Param.NAME);
		youngsEtAl1997AttenRel.getParameter(MagParam.NAME).setValue(6.0);
		youngsEtAl1997AttenRel.getParameter(DistanceRupParameter.NAME).setValue(10.0);
		youngsEtAl1997AttenRel.getParameter(Vs30_Param.NAME).setValue(800.0);
		youngsEtAl1997AttenRel.getParameter(TectonicRegionTypeParam.NAME).
			setValue(TectonicRegionType.SUBDUCTION_INTERFACE.toString());
		youngsEtAl1997AttenRel.getParameter(FocalDepthParam.NAME).setValue(20.0);
		
		youngsEtAl1997AttenRel.getParameter(PeriodParam.NAME).setValue(3.0);
		double mean3s = youngsEtAl1997AttenRel.getMean();
		
		youngsEtAl1997AttenRel.getParameter(PeriodParam.NAME).setValue(4.0);
		double mean4s = youngsEtAl1997AttenRel.getMean();
		
		assertTrue(mean4s  == mean3s / 0.399);
	}
	
	private void validateAgainstTable(final String tectonicRegionType,
			final double hypocentralDepth, final int periodIndex,
			final double mag, final String siteTypeParam,
			final int expectedResultIndex, final double[][] table) {
		for (int i = 0; i < table.length; i++) {

			double distance = table[i][0];
			double predicted = Math.exp(youngsEtAl1997AttenRel.
					getMean(periodIndex, mag, distance, tectonicRegionType,
							VS_30, hypocentralDepth));
			double expected = table[i][expectedResultIndex];
			double percentageDifference = Math.abs((expected - predicted)
					/ expected) * 100;

			String msg = "distance: " + distance + ", magnitude: " + mag
					+ ", site type param: " + siteTypeParam + ", tectonic region type: "
					+ tectonicRegionType + ", hypocentral depth: "
					+ hypocentralDepth + ", expected: " + expected
					+ ", predicted: " + predicted + ",percentage difference: "
					+ percentageDifference;
			assertTrue(msg, percentageDifference < TOLERANCE);
		}
	}

	private ArbitrarilyDiscretizedFunc setUpHazardCurve() {
		ArbitrarilyDiscretizedFunc hazCurve = new ArbitrarilyDiscretizedFunc();
		hazCurve.set(0.005, 0.0);
		hazCurve.set(0.007, 0.0);
		hazCurve.set(0.0098, 0.0);
		hazCurve.set(0.0137, 0.0);
		hazCurve.set(0.0192, 0.0);
		hazCurve.set(0.0269, 0.0);
		hazCurve.set(0.0376, 0.0);
		return hazCurve;
	}
	
	@Override
	public void parameterChangeWarning(ParameterChangeWarningEvent event) {
	}
	
}