package org.opensha.sha.imr.attenRelImpl;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;

import org.opensha.commons.data.NamedObjectAPI;
import org.opensha.commons.data.Site;
import org.opensha.commons.exceptions.ParameterException;
import org.opensha.commons.param.DoubleConstraint;
import org.opensha.commons.param.DoubleDiscreteConstraint;
import org.opensha.commons.param.StringConstraint;
import org.opensha.commons.param.event.ParameterChangeEvent;
import org.opensha.commons.param.event.ParameterChangeListener;
import org.opensha.commons.param.event.ParameterChangeWarningListener;
import org.opensha.sha.earthquake.EqkRupture;
import org.opensha.sha.imr.AttenuationRelationship;
import org.opensha.sha.imr.ScalarIntensityMeasureRelationshipAPI;
import org.opensha.sha.imr.attenRelImpl.constants.YoungsEtAl1997Constants;
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
 * <b>Title:</b> YoungsEtAl_1997_AttenRel
 * <p>
 * <b>Description:</b> Class implementing attenuation relationship described in:
 * "Strong Ground Motion Attenuation Relationships for Subduction Zone
 * Earthquakes", R.R.Youngs, S.J. Chiou, W.J. Silva, J. R. Humphrey,
 * Seismological Research Letters, Volume 68, Number 1, January/February 1997.
 * <p>
 * Supported Intensity-Measure Parameters:
 * <p>
 * <UL>
 * <LI>pgaParam - Peak Ground Acceleration
 * <LI>saParam - Response Spectral Acceleration (5% damping)
 * </UL>
 * <p>
 * Other Independent Parameters:
 * <p>
 * <UL>
 * <LI>magParam - moment magnitude (Mw)
 * <LI>distanceRupParam - closest distance to rupture surface (km)
 * <LI>vs30Param - shear wave velocity (m/s) averaged over the top 30 m of the
 * soil profile; The model assumes the following classification: vs30 > 760 ->
 * Rock, NEHRP A/B; vs30 <=760 -> soil, NEHRP C;
 * <LI>tectonicRegionTypeParam - interface or intraslab
 * <LI>focalDepthParam - depth to the earthquake rupture hypocenter (km)
 * <LI>componentParam - average horizontal
 * <LI>stdDevTypeParam - total, none
 * </UL>
 * <p>
 * 
 * <p>
 * 
 * Verification - This model has been validated (see
 * {@link YoungsEtAl_1997_test}) against tables provided by Céline Beauval
 * (<celine.beauval@obs.ujf-grenoble.fr>) using mathSHA. Tests were implemented
 * to check median PGA and SA (1Hz) for interface and intraslab events, at
 * different magnitude (8.8, 8.0, 7.0 for interface and 8.0, 7.5, 7.0 for
 * intraslab), for rock site type.
 * 
 * </p>
 * 
 ** 
 * @author L. Danciu, M.Pagani, D. Monelli
 * @version 1.0, December 2010
 */
public class YoungsEtAl_1997_AttenRel extends AttenuationRelationship implements
		ScalarIntensityMeasureRelationshipAPI, NamedObjectAPI,
		ParameterChangeListener {

	/** Short name. **/
	public static final String SHORT_NAME = "YoungsEtAl1997";

	/** Full name. **/
	public static final String NAME = "Youngs et. al. 1997";

	/** Version number. **/
	private static final long serialVersionUID = 1234567890987654353L;

	/** Moment magnitude (Mw). **/
	private double mag;

	/** Tectonic region type. **/
	private String tecRegType;

	/** Focal Depth. **/
	private double focalDepth;

	/** Vs30 **/
	private double vs30;

	/** Closest distance to rupture. **/
	private double rRup;

	/** Standard deviation type. **/
	private String stdDevType;

	/** Map period-value/period-index (for Rock). */
	private HashMap<Double, Integer> indexFromPerHashMapRock;

	/** Map period-value/period-index (for Soil). */
	private HashMap<Double, Integer> indexFromPerHashMapSoil;

	/** Period index. **/
	private int iper;

	/** For issuing warnings. **/
	private transient ParameterChangeWarningListener warningListener = null;

	/**
	 * Construct attenuation relationship. Initialize parameters and parameter
	 * lists.
	 */
	public YoungsEtAl_1997_AttenRel(
			ParameterChangeWarningListener warningListener) {

		// creates exceedProbParam
		super();

		this.warningListener = warningListener;

		initSupportedIntensityMeasureParams();

		// init period index - period value map for rock
		indexFromPerHashMapRock = new HashMap<Double, Integer>();
		for (int i = 1; i < YoungsEtAl1997Constants.PERIOD_ROCK.length; i++) {
			indexFromPerHashMapRock.put(new Double(
					YoungsEtAl1997Constants.PERIOD_ROCK[i]), new Integer(i));
		}
		// init period index - period value map for soil
		indexFromPerHashMapSoil = new HashMap<Double, Integer>();
		for (int i = 1; i < YoungsEtAl1997Constants.PERIOD_SOIL.length; i++) {
			indexFromPerHashMapSoil.put(new Double(
					YoungsEtAl1997Constants.PERIOD_SOIL[i]), new Integer(i));
		}

		initEqkRuptureParams();
		initSiteParams();
		initPropagationEffectParams();
		initOtherParams();
		initIndependentParamLists();
		initParameterEventListeners();
	}

	/**
	 * Creates the three supported IM parameters (PGA, PGV and SA), as well as
	 * the independenParameters of SA (periodParam and dampingParam) and adds
	 * them to the supportedIMParams list. Makes the parameters noneditable.
	 */
	protected void initSupportedIntensityMeasureParams() {

		// set supported periods for spectral acceleration
		// (supported periods are those for soil, same as rock except
		// the 4s period)
		DoubleDiscreteConstraint periodConstraint = new DoubleDiscreteConstraint();
		for (int i = 1; i < YoungsEtAl1997Constants.PERIOD_SOIL.length; i++) {
			periodConstraint.addDouble(new Double(
					YoungsEtAl1997Constants.PERIOD_SOIL[i]));
		}
		periodConstraint.setNonEditable();
		// set period param (default is 1s)
		saPeriodParam = new PeriodParam(periodConstraint);

		// set damping parameter. Empty constructor set damping
		// factor to 5 %
		saDampingParam = new DampingParam();

		// initialize spectral acceleration parameter (units: g)
		saParam = new SA_Param(saPeriodParam, saDampingParam);
		saParam.setNonEditable();

		// initialize peak ground acceleration parameter (units: g):
		pgaParam = new PGA_Param();
		pgaParam.setNonEditable();

		// add the warning listeners:
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
		magParam = new MagParam(YoungsEtAl1997Constants.MAG_WARN_MIN,
				YoungsEtAl1997Constants.MAG_WARN_MAX);

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
	 * Creates the Propagation Effect parameters and adds them to the
	 * propagationEffectParams list. Makes the parameters non-editable.
	 */
	protected void initPropagationEffectParams() {

		distanceRupParam = new DistanceRupParameter(0.0);
		distanceRupParam.addParameterChangeWarningListener(warningListener);
		DoubleConstraint warndistance = new DoubleConstraint(new Double(0.00),
				YoungsEtAl1997Constants.DISTANCE_RUP_WARN_MAX);
		warndistance.setNonEditable();
		distanceRupParam.setWarningConstraint(warndistance);
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
		stdDevTypeConstraint.setNonEditable();
		stdDevTypeParam = new StdDevTypeParam(stdDevTypeConstraint);

		// component Parameter
		StringConstraint constraint = new StringConstraint();
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
	public final void parameterChange(final ParameterChangeEvent e) {

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
	public final void resetParameterEventListeners() {
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
	public final void setEqkRupture(final EqkRupture eqkRupture) {

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
	public final void setSite(final Site site) {
		vs30Param.setValueIgnoreWarning((Double) site.getParameter(
				Vs30_Param.NAME).getValue());
		this.site = site;
		setPropagationEffectParams();
	}

	/**
	 * This sets the site and eqkRupture, and the related parameters, from the
	 * propEffect object passed in.
	 */

	public final void setPropagationEffectParams() {
		if ((this.site != null) && (this.eqkRupture != null)) {
			distanceRupParam.setValue(eqkRupture, site);
		}
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
	 * Set period index.
	 */
	protected void setPeriodIndex() throws ParameterException {

		if (im.getName().equalsIgnoreCase(PGA_Param.NAME)) {
			iper = 0;
		} else {
			double period = saPeriodParam.getValue();
			// rock
			if (vs30 > 760.0) {
				if (period > YoungsEtAl1997Constants.PERIOD_ROCK[YoungsEtAl1997Constants.PERIOD_ROCK.length - 1]) {
					throw new RuntimeException(period
							+ " not supported for vs30 > 760.0");
				} else {
					iper = ((Integer) indexFromPerHashMapRock.get(period))
							.intValue();
				}
			}
			// soil
			else {
				iper = ((Integer) indexFromPerHashMapSoil.get(period))
						.intValue();
			}
		}
	}

	/**
	 * Return mean.
	 */
	public double getMean() {

		if (rRup > USER_MAX_DISTANCE) {
			return VERY_SMALL_MEAN;
		}

		setPeriodIndex();

		return getMean(iper, mag, rRup, tecRegType, vs30, focalDepth);
	}

	/**
	 * Return standard deviation.
	 */
	public double getStdDev() {

		setPeriodIndex();

		double magnitude = mag;
		if (magnitude > 8.0) {
			magnitude = 8.0;
		}

		return getStdDev(iper, magnitude, stdDevType);
	}

	/**
	 * Compute mean.
	 */
	public double getMean(int iper, double mag, double rRup, String tecRegType,
			double vs30, double hypoDep) {

		double Zt, mean, lnY;
		if (tecRegType.equals(TectonicRegionType.SUBDUCTION_INTERFACE
				.toString())) {
			Zt = 0;
		} else {
			Zt = 1;
		}

		if (vs30 > 760) {
			// rock
			lnY = YoungsEtAl1997Constants.A1_ROCK
					+ YoungsEtAl1997Constants.A2_ROCK
					* mag
					+ YoungsEtAl1997Constants.C1_ROCK[iper]
					+ YoungsEtAl1997Constants.C2_ROCK[iper]
					* (Math.pow(YoungsEtAl1997Constants.A3_ROCK - mag, 3))
					+ YoungsEtAl1997Constants.C3_ROCK[iper]
					* Math.log(rRup + YoungsEtAl1997Constants.A4_ROCK
							* Math.exp(YoungsEtAl1997Constants.A5_ROCK * mag))
					+ YoungsEtAl1997Constants.A6_ROCK * hypoDep
					+ YoungsEtAl1997Constants.A7_ROCK * Zt;
			// correct for T = 4 s
			if (YoungsEtAl1997Constants.PERIOD_ROCK[iper]==4.00){
				lnY /= YoungsEtAl1997Constants.T3sec_TO_T4sec_DECAYFACTOR;
			}
		} else {
			// soil
			lnY = YoungsEtAl1997Constants.A1_SOIL
					+ YoungsEtAl1997Constants.A2_SOIL
					* mag
					+ YoungsEtAl1997Constants.C1_SOIL[iper]
					+ YoungsEtAl1997Constants.C2_SOIL[iper]
					* (Math.pow(YoungsEtAl1997Constants.A3_SOIL - mag, 3))
					+ YoungsEtAl1997Constants.C3_SOIL[iper]
					* Math.log(rRup + YoungsEtAl1997Constants.A4_SOIL
							* Math.exp(YoungsEtAl1997Constants.A5_SOIL * mag))
					+ YoungsEtAl1997Constants.A6_SOIL * hypoDep
					+ YoungsEtAl1997Constants.A7_SOIL * Zt;
		}

		return lnY;
	}

	/**
	 * Compute standard deviation.
	 */
	private double getStdDev(int iper, double mag, String stdDevType) {
		if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_NONE)) {
			return 0;
		}
		else if (stdDevType.equals(StdDevTypeParam.STD_DEV_TYPE_TOTAL)){
			double sigmaTotal = YoungsEtAl1997Constants.C4_SOIL[iper]
					+ YoungsEtAl1997Constants.C5_SOIL[iper] * mag;
			return (sigmaTotal);
		}
		throw new RuntimeException("Standard deviation type: "+stdDevType+" not supported");
	}

	/**
	 * Returns attenuation relationship full name.
	 */
	public String getName() {
		return NAME;
	}

	/**
	 * Returns attenuation relationship short name.
	 */
	public String getShortName() {
		return SHORT_NAME;
	}

	/**
	 * Provides URL with attenuation relationship info. The method currently
	 * returns null because no url has been created.
	 */
	public URL getInfoURL() throws MalformedURLException {
		return null;
	}
}
