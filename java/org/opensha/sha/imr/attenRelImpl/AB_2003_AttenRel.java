package org.opensha.sha.imr.attenRelImpl;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import org.opensha.commons.data.NamedObjectAPI;
import org.opensha.commons.data.Site;
import org.opensha.commons.exceptions.InvalidRangeException;
import org.opensha.commons.param.DoubleConstraint;
import org.opensha.commons.param.DoubleDiscreteConstraint;
import org.opensha.commons.param.StringConstraint;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.attenRelImpl.constants.AB2003Constants;
import org.opensha.sha.imr.param.EqkRuptureParams.FocalDepthParam;
import org.opensha.sha.imr.param.EqkRuptureParams.MagParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.DampingParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.PGA_Param;
import org.opensha.sha.imr.param.IntensityMeasureParams.PeriodParam;
import org.opensha.sha.imr.param.IntensityMeasureParams.SA_Param;
import org.opensha.sha.imr.param.OtherParams.ComponentParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncLevelParam;
import org.opensha.sha.imr.param.OtherParams.SigmaTruncTypeParam;
import org.opensha.sha.imr.param.OtherParams.StdDevTypeParam;
import org.opensha.sha.imr.param.OtherParams.TectonicRegionTypeParam;
import org.opensha.sha.imr.param.PropagationEffectParams.DistanceRupParameter;
import org.opensha.sha.imr.param.SiteParams.Vs30_Param;
import org.opensha.sha.util.TectonicRegionType;

/**
 * <b>Title:</b> AB_2003_AttenRel
 * <p>
 * <b>Description:</b> Class implementing attenuation relationship described in:
 * "Empirical Ground Motion Relations for Subduction-Zone Earthquakes and their
 * application to Cascadia and other regions", Gail M.Atkinson and David M.
 * Boore, Bulletin of the Seismological Society of America, Vol. 93, No. 4, pp
 * 1703-1729, 2003. "Erratum to 'Empirical Ground Motion Relations for
 * Subduction-Zone Earthquakes and their application to Cascadia and other
 * regions'", Gail M. Atkinson and David M. Boore, Vol. 98, No. 5, pp.2567-2569,
 * 2008. The class implements the global model but not the corrections for
 * Japan/Cascadia.
 * <p>
 * Supported Intensity-Measure Parameters:
 * <p>
 * <UL>
 * <LI>pgaParam - Peak Ground Acceleration
 * <LI>saParam - Response Spectral Acceleration
 * </UL>
 * <p>
 * Other Independent Parameters:
 * <p>
 * <UL>
 * <LI>magParam - moment magnitude
 * <LI>distanceRupParam - closest distance to rupture surface
 * <LI>vs30Param - shear wave velocity (m/s) averaged over the top 30 m of the
 * soil profile; The model assumes the following classification: vs30 > 760 ->
 * NEHRP B; 360 < vs30 <=760 -> NEHRP C; 180 <= vs30 <= 360 -> NEHRP D; vs30 <
 * 180 -> NEHRP E;
 * <LI>tectonicRegionTypeParam - interface or intra slab
 * <LI>focalDepthParam - depth to the earthquake rupture hypocenter
 * <LI>componentParam - random horizontal component, average horizontal (added
 * for SHARE projects requirement, assuming equivalence of random horizontal and
 * average horizontal)
 * <LI>stdDevTypeParam - total, inter-event, intra-event, none
 * </UL>
 * <p>
 * 
 * <p>
 * 
 * Verification - This model has been validated (see {@link AB_2003_test})
 * against tables provided by Céline Beauval
 * (<celine.beauval@obs.ujf-grenoble.fr>) using mathSHA. Tests were implemented
 * to check median PGA and SA (1Hz) for interface and intraslab events, at
 * different magnitude (7.0,8.0.8.8), for different site types (NERPH B, C, D).
 * 
 * 
 * </p>
 * 
 ** 
 * @author L. Danciu, D. Monelli
 * @version 1.0, December 2010
 */
public class AB_2003_AttenRel extends AttenuationRelationship implements
		ScalarIntensityMeasureRelationshipAPI, NamedObjectAPI,
		ParameterChangeListener {

	/** Short name. */
	public static String SHORT_NAME = "AB2003";

	/** Full name. */
	public static String NAME = "Atkinson & Boore 2003";

	/** Version number. */
	private static long serialVersionUID = 1234567890987654353L;

	/** Moment magnitude. */
	private double mag;

	/** Tectonic region type. */
	private String tecRegType;

	/** Focal depth. */
	private double focalDepth;

	/** Vs 30. */
	private double vs30;

	/** Closest distance to rupture. */
	private double rRup;

	/** Standard deviation type. */
	private String stdDevType;

	/** Map period-value/period-index. */
	private HashMap<Double, Integer> indexFromPerHashMap;

	/** Period index. */
	private int iper;

	/** For issuing warnings. */
	private transient ParameterChangeWarningListener warningListener = null;

	/**
	 * Construct attenuation relationship. Initialize parameters and parameter
	 * lists.
	 */
	public AB_2003_AttenRel(ParameterChangeWarningListener warningListener) {

		// creates exceedProbParam
		super();

		this.warningListener = warningListener;

		initSupportedIntensityMeasureParams();

		indexFromPerHashMap = new HashMap<Double, Integer>();
		for (int i = 1; i < AB2003Constants.PERIOD.length; i++) {
			indexFromPerHashMap.put(new Double(AB2003Constants.PERIOD[i]),
					new Integer(i));
		}

		initEqkRuptureParams();
		initSiteParams();
		initPropagationEffectParams();
		initOtherParams();
		initIndependentParamLists();
		initParameterEventListeners();
	}

	/**
	 * Creates the two supported IM parameters (PGA and SA), as well as the
	 * independenParameters of SA (periodParam and dampingParam) and adds them
	 * to the supportedIMParams list. Makes the parameters non-editable.
	 */
	protected void initSupportedIntensityMeasureParams() {

		// set supported periods for spectral acceleration
		DoubleDiscreteConstraint periodConstraint = new DoubleDiscreteConstraint();
		for (int i = 1; i < AB2003Constants.PERIOD.length; i++) {
			periodConstraint.addDouble(new Double(AB2003Constants.PERIOD[i]));
		}
		periodConstraint.setNonEditable();
		// set period param (default is 1s, which is provided by AB2003 GMPE)
		saPeriodParam = new PeriodParam(periodConstraint);

		// set damping parameter. Empty constructor set damping
		// factor to 5 % (which is the one provided by AB2003 GMPE)
		saDampingParam = new DampingParam();

		// initialize spectral acceleration parameter (units: g)
		saParam = new SA_Param(saPeriodParam, saDampingParam);
		saParam.setNonEditable();

		// initialize peak ground acceleration parameter (units: g)
		pgaParam = new PGA_Param();
		pgaParam.setNonEditable();

		// add the warning listeners
		saParam.addParameterChangeWarningListener(warningListener);
		pgaParam.addParameterChangeWarningListener(warningListener);

		// put parameters in the supportedIMParams list
		supportedIMParams.clear();
		supportedIMParams.addParameter(saParam);
		supportedIMParams.addParameter(pgaParam);
	}

	/**
	 * Initialize earthquake rupture parameter (moment magnitude, tectonic
	 * region type, focal depth) and add to eqkRuptureParams list. Makes the
	 * parameters non-editable.
	 */
	protected void initEqkRuptureParams() {

		// moment magnitude (default 5.5)
		magParam = new MagParam(AB2003Constants.MAG_WARN_MIN,
				AB2003Constants.MAG_WARN_MAX);

		// tectonic region type
		StringConstraint options = new StringConstraint();
		options.addString(TectonicRegionType.SUBDUCTION_INTERFACE.toString());
		options.addString(TectonicRegionType.SUBDUCTION_SLAB.toString());
		tectonicRegionTypeParam = new TectonicRegionTypeParam(options,
				TectonicRegionType.SUBDUCTION_INTERFACE.toString());

		// focal depth (default zero km)
		focalDepthParam = new FocalDepthParam();

		eqkRuptureParams.clear();
		eqkRuptureParams.addParameter(magParam);
		eqkRuptureParams.addParameter(tectonicRegionTypeParam);
		eqkRuptureParams.addParameter(focalDepthParam);
	}

	/**
	 * Initialize site parameters (vs30) and adds it to the siteParams list.
	 * Makes the parameters non-editable.
	 */
	protected void initSiteParams() {

		// vs30 parameters (constrains are not set, default value to 760 m/s)
		vs30Param = new Vs30_Param();

		siteParams.clear();
		siteParams.addParameter(vs30Param);
	}

	/**
	 * Initialize Propagation Effect parameters (closest distance to rupture)
	 * and adds them to the propagationEffectParams list. Makes the parameters
	 * non-editable.
	 */
	protected void initPropagationEffectParams() {

		distanceRupParam = new DistanceRupParameter(
				AB2003Constants.DISTANCE_RUP_WARN_MIN);
		distanceRupParam.addParameterChangeWarningListener(warningListener);
		DoubleConstraint warn = new DoubleConstraint(new Double(0.00),
				AB2003Constants.DISTANCE_RUP_WARN_MAX);
		warn.setNonEditable();
		distanceRupParam.setWarningConstraint(warn);
		distanceRupParam.setNonEditable();

		propagationEffectParams.addParameter(distanceRupParam);
	}

	/**
	 * Initialize other Parameters (standard deviation type, component, sigma
	 * truncation type, sigma truncation level).
	 */
	protected void initOtherParams() {

		sigmaTruncTypeParam = new SigmaTruncTypeParam();
		sigmaTruncLevelParam = new SigmaTruncLevelParam();

		// stdDevType Parameter
		StringConstraint stdDevTypeConstraint = new StringConstraint();
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_TOTAL);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_NONE);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_INTER);
		stdDevTypeConstraint.addString(StdDevTypeParam.STD_DEV_TYPE_INTRA);
		stdDevTypeConstraint.setNonEditable();
		stdDevTypeParam = new StdDevTypeParam(stdDevTypeConstraint);

		// component Parameter
		StringConstraint constraint = new StringConstraint();
		constraint.addString(ComponentParam.COMPONENT_RANDOM_HORZ);
		constraint.addString(ComponentParam.COMPONENT_AVE_HORZ);
		constraint.setNonEditable();
		componentParam = new ComponentParam(constraint,
				ComponentParam.COMPONENT_AVE_HORZ);

		// add these to the list
		otherParams.clear();
		otherParams.addParameter(sigmaTruncTypeParam);
		otherParams.addParameter(sigmaTruncLevelParam);
		otherParams.addParameter(stdDevTypeParam);
		otherParams.addParameter(componentParam);
	}

	/**
	 * This creates the lists of independent parameters that the various
	 * dependent parameters (mean, standard deviation, exceedance probability,
	 * and IML at exceedance probability) depend upon. NOTE: these lists do not
	 * include anything about the intensity-measure parameters or any of their
	 * internal independentParamaters.
	 */
	protected void initIndependentParamLists() {

		// params that the mean depends upon
		meanIndependentParams.clear();
		meanIndependentParams.addParameter(magParam);
		meanIndependentParams.addParameter(tectonicRegionTypeParam);
		meanIndependentParams.addParameter(focalDepthParam);
		meanIndependentParams.addParameter(vs30Param);
		meanIndependentParams.addParameter(distanceRupParam);

		// params that the stdDev depends upon
		stdDevIndependentParams.clear();
		stdDevIndependentParams.addParameter(tectonicRegionTypeParam);
		stdDevIndependentParams.addParameter(stdDevTypeParam);

		// params that the exceed. prob. depends upon
		exceedProbIndependentParams.clear();
		exceedProbIndependentParams.addParameterList(meanIndependentParams);
		exceedProbIndependentParams.addParameter(stdDevTypeParam);
		exceedProbIndependentParams.addParameter(sigmaTruncTypeParam);
		exceedProbIndependentParams.addParameter(sigmaTruncLevelParam);

		// params that the IML at exceed. prob. depends upon
		imlAtExceedProbIndependentParams
				.addParameterList(exceedProbIndependentParams);
		imlAtExceedProbIndependentParams.addParameter(exceedProbParam);
	}

	/**
	 * Adds the parameter change listeners. This allows to listen to when-ever
	 * the parameter is changed.
	 */
	protected void initParameterEventListeners() {

		// earthquake rupture params
		magParam.addParameterChangeListener(this);
		tectonicRegionTypeParam.addParameterChangeListener(this);
		focalDepthParam.addParameterChangeListener(this);

		// site params
		vs30Param.addParameterChangeListener(this);

		// propagation effect param
		distanceRupParam.addParameterChangeListener(this);

		// standard deviation type param
		stdDevTypeParam.addParameterChangeListener(this);
	}

	/**
	 * This listens for parameter changes and updates the primitive parameters
	 * accordingly
	 */
	public void parameterChange(ParameterChangeEvent e) {

		String pName = e.getParameterName();
		Object val = e.getNewValue();

		if (pName.equals(MagParam.NAME)) {
			mag = ((Double) val).doubleValue();
		} else if (pName.equals(TectonicRegionTypeParam.NAME)) {
			tecRegType = (String) val;
		} else if (pName.equals(FocalDepthParam.NAME)) {
			focalDepth = ((Double) val).doubleValue();
		} else if (pName.equals(Vs30_Param.NAME)) {
			vs30 = ((Double) val).doubleValue();
		} else if (pName.equals(DistanceRupParameter.NAME)) {
			rRup = ((Double) val).doubleValue();
		} else if (pName.equals(StdDevTypeParam.NAME)) {
			stdDevType = (String) val;
		}
	}

	/**
	 * Allows to reset the change listeners on the parameters.
	 */
	public void resetParameterEventListeners() {
		magParam.removeParameterChangeListener(this);
		tectonicRegionTypeParam.removeParameterChangeListener(this);
		focalDepthParam.removeParameterChangeListener(this);
		vs30Param.removeParameterChangeListener(this);
		distanceRupParam.removeParameterChangeListener(this);
		stdDevTypeParam.removeParameterChangeListener(this);
		this.initParameterEventListeners();
	}

	/**
	 * This sets the eqkRupture related parameters (moment magnitude, tectonic
	 * region type, focal depth) based on the eqkRupture passed in. The
	 * internally held eqkRupture object is also set as that passed in. Warning
	 * constrains on magnitude and focal depth are ignored.
	 */
	public void setEqkRupture(EqkRupture eqkRupture)
			throws InvalidRangeException {
		magParam.setValueIgnoreWarning(new Double(eqkRupture.getMag()));

		if (eqkRupture.getTectRegType() != null) {
			tectonicRegionTypeParam.setValue(eqkRupture.getTectRegType()
					.toString());
		} else {
			throw new RuntimeException("Tectonic region type not set in "
					+ " earthquake rupture");
		}

		if (eqkRupture.getHypocenterLocation() != null) {
			focalDepthParam.setValueIgnoreWarning(new Double(eqkRupture
					.getHypocenterLocation().getDepth()));
		} else {
			throw new RuntimeException("Hypocenter location not set in"
					+ " earthquake rupture");
		}

		this.eqkRupture = eqkRupture;

		setPropagationEffectParams();
	}

	/**
	 * This sets the site-related parameter (vs30) based on what is in the Site
	 * object passed in. This also sets the internally held Site object as that
	 * passed in.
	 */
	public void setSite(Site site) {
		vs30Param.setValueIgnoreWarning((Double) site.getParameter(
				Vs30_Param.NAME).getValue());
		this.site = site;
		setPropagationEffectParams();
	}

	/**
	 * This sets the site and eqkRupture, and the related parameters, from the
	 * propEffect object passed in.
	 */
	public void setPropagationEffectParams() {
		if ((this.site != null) && (this.eqkRupture != null)) {
			distanceRupParam.setValue(eqkRupture, site);
		}
	}

	/**
	 * Set period index.
	 */
	protected void setPeriodIndex() {
		if (im.getName().equalsIgnoreCase(PGA_Param.NAME)) {
			iper = 0;
		} else {
			iper = ((Integer) indexFromPerHashMap.get(saPeriodParam.getValue()))
					.intValue();
		}
	}

	/**
	 * Compute mean. Applies correction for periods = 2.5 and 5 Hz (for
	 * interface) as for Atkinson and Boore 2008 Erratum.
	 */
	public double getMean() {

		if (rRup > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		}

		double mean = Double.NaN;

		// if period corresponds to 5 Hz (PERIOD[3]) or 2.5 (PERIOD[4]) Hz
		// and tectonic region type is interface applies correction.
		if (im.getName().equalsIgnoreCase(SA_Param.NAME) 
				&& (saPeriodParam.getValue() == AB2003Constants.PERIOD[3] || saPeriodParam
				.getValue() == AB2003Constants.PERIOD[4])
				&& tecRegType
						.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_INTERFACE
								.toString())) {
			// compute log SA for 2.5 and 5 Hz and return a weighted sum
			int indexPeriod5Hz = ((Integer) indexFromPerHashMap
					.get(AB2003Constants.PERIOD[3])).intValue();
			int indexPeriod25Hz = ((Integer) indexFromPerHashMap
					.get(AB2003Constants.PERIOD[4])).intValue();
			double logSA5Hz = getMean(indexPeriod5Hz, mag, rRup, vs30,
					tecRegType, focalDepth);
			double logSA25Hz = getMean(indexPeriod25Hz, mag, rRup, vs30,
					tecRegType, focalDepth);
			if (saPeriodParam.getValue() == AB2003Constants.PERIOD[3]) {
				mean = AB2003Constants.CORRECTION_WEIGHTS[0] * logSA5Hz
						+ AB2003Constants.CORRECTION_WEIGHTS[1] * logSA25Hz;
				return mean;
			} else if (saPeriodParam.getValue() == AB2003Constants.PERIOD[4]) {
				mean = AB2003Constants.CORRECTION_WEIGHTS[0] * logSA25Hz
						+ AB2003Constants.CORRECTION_WEIGHTS[1] * logSA5Hz;
				return mean;
			}
		} else {
			setPeriodIndex();
			return getMean(iper, mag, rRup, vs30, tecRegType, focalDepth);
		}
		return mean;
	}

	/**
	 * Compute standard deviation.
	 */
	public double getStdDev() {

		setPeriodIndex();

		return getStdDev(iper, stdDevType, tecRegType);
	}

	/**
	 * Allows the user to set the default parameter values for the selected
	 * Attenuation Relationship.
	 */
	public void setParamDefaults() {

		magParam.setValueAsDefault();
		tectonicRegionTypeParam.setValueAsDefault();
		focalDepthParam.setValueAsDefault();
		vs30Param.setValueAsDefault();
		distanceRupParam.setValueAsDefault();
		saPeriodParam.setValueAsDefault();
		saDampingParam.setValueAsDefault();
		saParam.setValueAsDefault();
		pgaParam.setValueAsDefault();
		stdDevTypeParam.setValueAsDefault();
		sigmaTruncTypeParam.setValueAsDefault();
		sigmaTruncLevelParam.setValueAsDefault();
		componentParam.setValueAsDefault();
	}

	/**
	 * Get the name of this IMR.
	 */
	public String getName() {
		return NAME;
	}

	/**
	 * Returns the Short Name of each AttenuationRelationship
	 * 
	 */
	public String getShortName() {
		return SHORT_NAME;
	}

	/**
	 * Compute mean (natural logarithm of median ground motion).
	 */
	public double getMean(int iper, double mag, double rRup, double vs30,
			String tecRegType, double hypoDep) {

		hypoDep = capHypocentralDepth(hypoDep);

		mag = capMagnitude(mag, tecRegType);

		double rockResponse = computeRockResponse(tecRegType, iper, mag,
				hypoDep, rRup);

		double soilResponse = computeSoilResponse(tecRegType, iper, mag,
				hypoDep, rRup, vs30);

		double logY = rockResponse + soilResponse;
		logY *= AB2003Constants.LOG10_2_LN;

		double mean = Math.log(Math.exp(logY)
				* AB2003Constants.CMS_TO_G_CONVERSION_FACTOR);
		
		// apply correction for T = 4 s
		if(AB2003Constants.PERIOD[iper] == 4.0){
			mean /= AB2003Constants.T3sec_TO_T4sec_DECAYFACTOR;
		}
		
		return mean;
	}

	private double capHypocentralDepth(double hypoDep) {
		if (hypoDep > AB2003Constants.THRESHOLD_HYPO_DEPTH) {
			hypoDep = AB2003Constants.THRESHOLD_HYPO_DEPTH;
		}
		return hypoDep;
	}

	private double capMagnitude(double mag, String tecRegType) {
		double thresholdMag = Double.NaN;
		if (tecRegType.equals(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			thresholdMag = AB2003Constants.THRESHOLD_MAG_INTERFACE;
		} else if (tecRegType.equals(TectonicRegionType.SUBDUCTION_SLAB
				.toString())) {
			thresholdMag = AB2003Constants.THRESHOLD_MAG_INTRASLAB;
		}
		if (mag >= thresholdMag) {
			mag = thresholdMag;
		}
		return mag;
	}

	private double computeRockResponse(String tecRegType, int periodIndex,
			double mag, double hypoDep, double rRup) {
		double delta = AB2003Constants.NEAR_SOURCE_SATURATION_FACTOR1
				* Math.pow(10, AB2003Constants.NEAR_SOURCE_SATURATION_FACTOR2
						* mag);
		double R = Math.sqrt(rRup * rRup + delta * delta);
		double g = computedGeometricSpreadingFactor(tecRegType, mag);
		double rockResponse = Double.NaN;
		if (tecRegType.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			rockResponse = AB2003Constants.INTER_C1[periodIndex]
					+ AB2003Constants.INTER_C2[periodIndex] * mag
					+ AB2003Constants.INTER_C3[periodIndex] * hypoDep
					+ AB2003Constants.INTER_C4[periodIndex] * R - g
					* Math.log10(R);
		} else if (tecRegType
				.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_SLAB.toString())) {
			rockResponse = AB2003Constants.INTRA_C1[periodIndex]
					+ AB2003Constants.INTRA_C2[periodIndex] * mag
					+ AB2003Constants.INTRA_C3[periodIndex] * hypoDep
					+ AB2003Constants.INTRA_C4[periodIndex] * R - g
					* Math.log10(R);
		}
		return rockResponse;
	}

	/**
	 * Compute geometric spreading factor (i.e. 10^(f1-f2*mag))
	 */
	private double computedGeometricSpreadingFactor(String tecRegType,
			double mag) {
		double f1 = Double.NaN;
		double f2 = Double.NaN;
		if (tecRegType.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			f1 = AB2003Constants.INTERFACE_GEOM_SPREAD_FACTOR1;
			f2 = AB2003Constants.INTERFACE_GEOM_SPREAD_FACTOR2;
		} else if (tecRegType
				.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_SLAB.toString())) {
			f1 = AB2003Constants.INTRASLAB_GEOM_SPREAD_FACTOR1;
			f2 = AB2003Constants.INTRASLAB_GEOM_SPREAD_FACTOR2;
		}
		return Math.pow(10, f1 - f2 * mag);
	}

	private double computeSoilResponse(String tecRegType, int periodIndex,
			double mag, double hypoDep, double rRup, double vs30) {
		int indexPga = 0;
		double pgaOnRock = Math.pow(10, computeRockResponse(tecRegType,
				indexPga, mag, hypoDep, rRup));
		double sl = computeSoilLinearityTerm(periodIndex, pgaOnRock);
		double[] s = computeSiteTermCorrection(vs30);
		double soilResponse = Double.NaN;
		if (tecRegType.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			soilResponse = s[0] * AB2003Constants.INTER_C5[periodIndex] * sl
					+ s[1] * AB2003Constants.INTER_C6[periodIndex] * sl + s[2]
					* AB2003Constants.INTER_C7[periodIndex] * sl;
		} else if (tecRegType
				.equalsIgnoreCase(TectonicRegionType.SUBDUCTION_SLAB.toString())) {
			soilResponse = s[0] * AB2003Constants.INTRA_C5[periodIndex] * sl
					+ s[1] * AB2003Constants.INTRA_C6[periodIndex] * sl + s[2]
					* AB2003Constants.INTRA_C7[periodIndex] * sl;
		}
		return soilResponse;
	}

	private double computeSoilLinearityTerm(int iper, double PGArx) {
		double sl = Double.NaN;
		if (PGArx <= 100.0 || AB2003Constants.FREQ[iper] <= 1.0) {
			sl = 1.0;
		} else if (PGArx > 100.0 && PGArx < 500.0) {
			if (AB2003Constants.FREQ[iper] > 1.0
					&& AB2003Constants.FREQ[iper] < 2.0) {
				sl = 1.0 - (AB2003Constants.FREQ[iper] - 1.0) * (PGArx - 100.0)
						/ 400.0;
			} else if (AB2003Constants.FREQ[iper] > 2.0) {
				sl = 1.0 - (PGArx - 100.0) / 400.0;
			}
		} else if (PGArx >= 500.0) {
			if (AB2003Constants.FREQ[iper] > 1.0
					&& AB2003Constants.FREQ[iper] < 2.0) {
				sl = 1.0 - (AB2003Constants.FREQ[iper] - 1.0);
			} else if (AB2003Constants.FREQ[iper] > 2.0) {
				sl = 0.0;
			}
		}
		return sl;
	}

	private double[] computeSiteTermCorrection(double vs30) {
		double[] s = new double[3];
		if (vs30 > AB2003Constants.NEHRP_C_UPPER_BOUND) {
			s[0] = 0.0;
			s[1] = 0.0;
			s[2] = 0.0;
		} else if (vs30 > AB2003Constants.NEHRP_D_UPPER_BOUND
				&& vs30 <= AB2003Constants.NEHRP_C_UPPER_BOUND) {
			s[0] = 1.0;
			s[1] = 0.0;
			s[2] = 0.0;
		} else if (vs30 >= AB2003Constants.NEHRP_E_UPPER_BOUND
				&& vs30 <= AB2003Constants.NEHRP_D_UPPER_BOUND) {
			s[0] = 0.0;
			s[1] = 1.0;
			s[2] = 0.0;
		} else if (vs30 < AB2003Constants.NEHRP_E_UPPER_BOUND) {
			s[0] = 0.0;
			s[1] = 0.0;
			s[2] = 1.0;
		}
		return s;
	}

	/**
	 * This gets the standard deviation for specific parameter settings.
	 */
	public double getStdDev(int iper, String stdDevType, String tecRegType) {

		if (tecRegType.equals(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_TOTAL)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTER_TOTAL_STD[iper];
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_NONE)) {
				return 0;
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTER)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTER_INTEREVENT_STD[iper];
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTRA)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTER_INTRAEVENT_STD[iper];
			} else {
				throw new RuntimeException("Standard deviation type not recognized");
			}
		} else if (tecRegType.equals(TectonicRegionType.SUBDUCTION_SLAB
				.toString())) {
			if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_TOTAL)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTRA_TOTAL_STD[iper];
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_NONE)) {
				return 0;
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTER)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTRA_INTEREVENT_STD[iper];
			} else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_INTRA)) {
				return AB2003Constants.LOG10_2_LN
						* AB2003Constants.INTRA_INTRAEVENT_STD[iper];
			} else {
				throw new RuntimeException("Standard deviation type not recognized");
			}
		}
		throw new RuntimeException("Tectonic region type not recognized");
	}

	/**
	 * This provides a URL where more info on this model can be obtained. It
	 * currently returns null because no URL has been set up.
	 */
	public URL getInfoURL() throws MalformedURLException {
		return null;
	}

}
