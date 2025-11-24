/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SmoothIMK.cpp,v $

// Written: S. A. Jalali 10/2019
// Adding Cyclic and in-cycle deterioration modes to steel02 UniaxialMaterial

#include <math.h>

#include <stdlib.h>
#include <SmoothIMK.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <MaterialResponse.h>

#include <elementAPI.h>
#include <OPS_Globals.h>
#include "classTags.h"

static int numThisCall = 0;
void printSyntax()
{
	opserr << "------ SmoothIMK unaxialMaterial -------\n";
	opserr << "-------Syntax:\n";
	opserr << "-------UniaxialMaterial SmoothIMK $matTag,\n";
	opserr << "				 -posBackBone pd1 pf1 pd2 pf2 pd3 pf3 ... <-gap gap=0>\n";
	opserr << "        <-negBackBone nd1 nf1 nd2 nf2 nd3 nf3 ... <-gap gap=0>>\n";
	opserr << "        <-sigInit sigInit>\n";
	opserr << "        <-deterioration gamaS cS gamaK cK>\n";
	opserr << "        <-transition r0 r1=26 r2=0.2> \n";
	opserr << "        <-posCyclic <-bilinEndAmp value=1>\n";
	opserr << "                   <-peakOriented>\n";
	opserr << "                   <-pinched alphaPinch pinchY <-nonlin betaPinch=0 epsPinch=1>>\n";
	opserr << "                   <-unloading unloadingStiffFac <stressPenetFac=0>>\n";
	opserr << "        >\n";
	opserr << "        <-negCyclic <-bilinEndAmp value=1>\n";
	opserr << "                   <-peakOriented>\n";
	opserr << "                   <-pinched alphaPinch pinchY <-nonlin betaPinch=0 epsPinch=1>>\n";
	opserr << "                   <-unloading unloadingStiffFac <stressPenetFac=0>>\n";
	opserr << "        >\n";
}
//bool findStr(std::vector<const char*> vec, const char* str)
//{
//	for (auto iter = vec.begin(); iter != vec.end(); iter++)
//		if (strcmp(*iter, str) == 0)
//			return true;
//	return false;
//}
void*
OPS_SmoothIMK()
{
	if (numThisCall == 0) {
		numThisCall = 1;
	}
	// Pointer to a uniaxial material that will be returned
	UniaxialMaterial* theMaterial = 0;
	// List of recognized option strings for SmoothIMK
	//const std::vector<const char*> smoothIMKOptions = {
	//		"-negBackBone",
	//		"-sigInit",
	//		"-gap",
	//		"-deterioration",
	//		"-transition",
	//		"-peakOriented",
	//		"-pinched",
	//		"-unloadingStiffFactor"
	//};

	int numData = 1;
	int tag = 0;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial SmoothIMK tag" << endln;
		printSyntax();
		return 0;
	}

	if (OPS_GetNumRemainingInputArgs() < 7) {
		opserr << "Invalid SmoothIMK #args for: " << tag << " see the syntax" << endln;
		printSyntax();
		return 0;
	}
	std::vector<double> pf, pd, nf, nd;
	double bilinEndAmpPos = 1.2, bilinEndAmpNeg = 1.2;
	double unloadingStiffFacPos = 1, unloadingStiffFacNeg = 1; //default values for unloading stiffness factors
	//default parameters:
	double gamaS = 1e6, cS = 1, gamaUE = 1e6, cUE = 1, r0 = 3, r1 = 26, r2 = 0.2, sigInit = 0;
	double alphaPinchPos = 0, betaPinchPos = 0, epsPinchPos = 1, pinchYPos = 1,
		alphaPinchNeg = 0, pinchYNeg = 1, betaPinchNeg = 0, epsPinchNeg = 1;
	double sigPenetFacP = 0, sigPenetFacN = 0;
	double gapP = 0, gapN = 0;
	int ruleP = 1, ruleN = 1;
	const char* option = "";
	double val;
	int n = 0;
	while (OPS_GetNumRemainingInputArgs() > 0)
	{
		option = OPS_GetString();
		if (strcmp(option, "-posBackBone") == 0) {
			n = 0;
			while (OPS_GetDoubleInput(&numData, &val) == 0)
			{
				n++;
				if (n % 2 == 1)
					pd.push_back(val);
				else
					pf.push_back(val);
			}
			n = pd.size();
			if (pf.size() != n)
			{
				opserr << "SmoothIMK:: Error: the number of positive force and displacement inputs are not the same for material with tag: " << tag << endln;
				return 0;
			}
			if (n < 2)
			{
				opserr << "SmoothIMK:: Error: at least two positive force and displacement inputs are required for material with tag: " << tag << endln;
				return 0;
			}
			for (int i = 0; i < n - 1; i++)
			{
				if (pd[i] >= pd[i + 1])
				{
					char buf[16];
					snprintf(buf, sizeof(buf), "pd%d >= pd%d", i, i + 1);
					opserr << "SmoothIMK:: Error: " << buf << "in material with tag: " << tag << endln;
					printSyntax();
					return 0;
				}
			}
			option = OPS_GetString();
			if (strcmp(option, "-gap") == 0)
			{
				if (OPS_GetDoubleInput(&numData, &gapP) != 0) {
					opserr << "SmoothIMK:: invalid positive gap for material : " << tag << endln;
					return 0;
				}
				gapN = -gapP;
			}
			else if (strcmp(option, "Invalid String Input!") != 0)
				OPS_ResetCurrentInputArg(-1);
		}
		else if (strcmp(option, "-negBackBone") == 0) {
			n = 0;
			while (OPS_GetDoubleInput(&numData, &val) == 0)
			{
				n++;
				if (n % 2 == 1)
					nd.push_back(val < 0 ? val : -val);
				else
					nf.push_back(val < 0 ? val : -val);
			}
			n = nd.size();
			if (nf.size() != n)
			{
				opserr << "SmoothIMK:: Error: the number of negative force and displacement inputs are not the same for material with tag: " << tag << endln;
				return 0;
			}
			if (n < 2)
			{
				opserr << "SmoothIMK:: Error: at least two negative force and displacement inputs are required for material with tag: " << tag << endln;
				return 0;
			}
			for (int i = 0; i < n - 1; i++)
			{
				if (nd[i] <= nd[i + 1])
				{
					char buf[16];
					snprintf(buf, sizeof(buf), "nd%d <= nd%d", i, i + 1);
					opserr << "SmoothIMK:: Error: " << buf << "in material with tag: " << tag << endln;
					printSyntax();
					return 0;
				}
			}
			option = OPS_GetString();
			if (strcmp(option, "-gap") == 0)
			{
				if (OPS_GetDoubleInput(&numData, &gapN) != 0) {
					opserr << "SmoothIMK:: invalid negative gap for material : " << tag << endln;
					return 0;
				}
				gapN = -fabs(gapN);
			}
			else if (strcmp(option, "Invalid String Input!") != 0)
				OPS_ResetCurrentInputArg(-1);
		}
		else if (strcmp(option, "-sigInit") == 0) {
			if (OPS_GetDoubleInput(&numData, &sigInit) != 0) {
				opserr << "SmoothIMK:: invalid sigInit for material : " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-deterioration") == 0) {
			if (OPS_GetDoubleInput(&numData, &gamaS) != 0) {
				opserr << "SmoothIMK:: invalid gamaS for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &cS) != 0) {
				opserr << "SmoothIMK:: invalid cS for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &gamaUE) != 0) {
				opserr << "SmoothIMK:: invalid gamaUloadE for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &cUE) != 0) {
				opserr << "SmoothIMK:: invalid cUloadE for material : " << tag << endln;
				return 0;
			}
		}
		else if (strcmp(option, "-transition") == 0) {
			if (OPS_GetDoubleInput(&numData, &r0) != 0) {
				opserr << "SmoothIMK:: invalid r0 for material : " << tag << endln;
				return 0;
			}
			if (OPS_GetDoubleInput(&numData, &r1) != 0) {
				//OPS_ResetCurrentInputArg(-1);
				r1 = 0, r2 = 0;
			}
			else {
				if (OPS_GetDoubleInput(&numData, &r2) != 0) {
					opserr << "SmoothIMK:: invalid r2 for material : " << tag << endln;
					return 0;
				}
			}
		}
		else if (strcmp(option, "-posCyclic") == 0)
		{
			while (OPS_GetNumRemainingInputArgs() > 0)
			{
				option = OPS_GetString();
				if (strcmp(option, "-peakOriented") == 0)
				{
					ruleP = 2;
					alphaPinchPos = 0, betaPinchPos = 0, epsPinchPos = 1, pinchYPos = 0;
				}
				else if (strcmp(option, "-bilinEndAmp") == 0)
				{
					if (OPS_GetDoubleInput(&numData, &bilinEndAmpPos) != 0)
					{
						opserr << "SmoothIMK:: invalid positive startAmp for peakOriented rule for material : " << tag << endln;
						return 0;
					}
				}
				else if (strcmp(option, "-unloading") == 0)
				{
					if (OPS_GetDoubleInput(&numData, &unloadingStiffFacPos) != 0)
					{
						opserr << "SmoothIMK:: invalid positive unloading stiffness factor for material: " << tag << endln;
						return 0;
					}
					OPS_GetDoubleInput(&numData, &sigPenetFacP);
				}
				else if (strcmp(option, "-pinched") == 0)
				{
					ruleP = 3;
					if (OPS_GetDoubleInput(&numData, &alphaPinchPos) != 0)
					{
						opserr << "SmoothIMK:: invalid positive pinching alpha for material : " << tag << endln;
						return 0;
					}
					if (OPS_GetDoubleInput(&numData, &pinchYPos) != 0)
					{
						opserr << "SmoothIMK:: invalid positive stress pinching factor for material : " << tag << endln;
						return 0;
					}
					option = OPS_GetString();
					betaPinchPos = 0, epsPinchPos = 1;
					if (strcmp(option, "-nonlin") == 0)
					{
						if (OPS_GetDoubleInput(&numData, &betaPinchPos) != 0)
						{
							opserr << "SmoothIMK:: invalid positive pinching factor power for material : " << tag << endln;
							return 0;
						}
						if (OPS_GetDoubleInput(&numData, &epsPinchPos) != 0)
						{
							opserr << "SmoothIMK:: invalid positive pinching power normalizer for material : " << tag << endln;
							return 0;
						}
					}
					else if (strcmp(option, "Invalid String Input!") != 0)
						OPS_ResetCurrentInputArg(-1);
				}
				else
				{
					if (strcmp(option, "Invalid String Input!") != 0)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
			}
		}
		else if (strcmp(option, "-negCyclic") == 0)
		{
			while (OPS_GetNumRemainingInputArgs() > 0)
			{
				option = OPS_GetString();
				if (strcmp(option, "-peakOriented") == 0)
				{
					ruleN = 2;
					alphaPinchNeg = 0, betaPinchNeg = 0, epsPinchNeg = 1, pinchYNeg = 0;
				}
				else if (strcmp(option, "-bilinEndAmp") == 0)
				{
					if (OPS_GetDoubleInput(&numData, &bilinEndAmpNeg) != 0)
					{
						opserr << "SmoothIMK:: invalid negative startAmp for peakOriented rule for material : " << tag << endln;
						return 0;
					}
				}
				else if (strcmp(option, "-unloading") == 0)
				{
					if (OPS_GetDoubleInput(&numData, &unloadingStiffFacNeg) != 0)
					{
						opserr << "SmoothIMK:: invalid negative unloading stiffness factor for material: " << tag << endln;
						return 0;
					}
					OPS_GetDoubleInput(&numData, &sigPenetFacN);
				}
				else if (strcmp(option, "-pinched") == 0)
				{
					ruleN = 3;
					if (OPS_GetDoubleInput(&numData, &alphaPinchNeg) != 0)
					{
						opserr << "SmoothIMK:: invalid negative pinching alpha for material : " << tag << endln;
						return 0;
					}
					if (OPS_GetDoubleInput(&numData, &pinchYNeg) != 0)
					{
						opserr << "SmoothIMK:: invalid negative stress pinching factor for material : " << tag << endln;
						return 0;
					}
					option = OPS_GetString();
					betaPinchNeg = 0, epsPinchNeg = 1;
					if (strcmp(option, "-nonlin") == 0)
					{
						if (OPS_GetDoubleInput(&numData, &betaPinchNeg) != 0)
						{
							opserr << "SmoothIMK:: invalid negative pinching factor power for material : " << tag << endln;
							return 0;
						}
						if (OPS_GetDoubleInput(&numData, &epsPinchNeg) != 0)
						{
							opserr << "SmoothIMK:: invalid negative pinching power normalizer for material : " << tag << endln;
							return 0;
						}
					}
					else if (strcmp(option, "Invalid String Input!") != 0)
						OPS_ResetCurrentInputArg(-1);
				}
				else
				{
					if (strcmp(option, "Invalid String Input!") != 0)
						OPS_ResetCurrentInputArg(-1);
					break;
				}
			}
		}
		else {
			if (strcmp(option, "Invalid String Input!") == 0)
			{
				//OPS_ResetCurrentInputArg(-1);
				char data[128];
				OPS_GetStringFromAll(data, 128);
				opserr << "SmoothIMK::Error: expected string switch but got non-string input: " << data << " for material with tag: " << tag << endln;
			}
			else
				opserr << "SmoothIMK::Error: Invalid option: " << option << endln;
			return 0;
		}
	}
	if (nd.size() == 0)
	{
		// No negative backbone was given, so we use the positive backbone for both directions
		nd.resize(pd.size());
		for (size_t i = 0; i < pd.size(); ++i) {
			nd[i] = -pd[i];
		}

		nf.resize(pf.size());
		for (size_t i = 0; i < pf.size(); ++i) {
			nf[i] = -pf[i];
		}
	}
	// Parsing was successful, allocate the material
	theMaterial = new SmoothIMK(tag, pd, pf, nd, nf,
		gamaS, cS, gamaUE, cUE, r0, r1, r2,
		ruleP, gapP, alphaPinchPos, pinchYPos, betaPinchPos, epsPinchPos, sigPenetFacP, bilinEndAmpPos, unloadingStiffFacPos,
		ruleN, gapN, alphaPinchNeg, pinchYNeg, betaPinchNeg, epsPinchNeg, sigPenetFacN, bilinEndAmpNeg, unloadingStiffFacNeg, sigInit);


	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type SmoothIMK\n";
		return 0;
	}

	return theMaterial;
}

double pinchDx(double x, double alphaPinch, double epsPinch, double betaPinch, double x0)
{
	if (betaPinch == 0)
		return alphaPinch * x;
	if (fabs(x) < 1e-6)
	{
		//avoid division by zero
		x = x > 0 ? 1e-6 : -1e-6;
	}
	return pow(alphaPinch, pow(epsPinch * x0 / x, betaPinch)) * x;
}

SmoothIMK::SmoothIMK(int tag,
	std::vector<double> _pd, std::vector<double> _pf,
	std::vector<double> _nd, std::vector<double> _nf,
	double _gamaS, double _cS,
	double _gamaUE, double _cUE,
	double _r0, double _r1, double _r2,
	int _cyclicRuleP, double _gapP, double _alphaPinchPos, double _pinchYPos,
	double _betaPinchPos, double _epsPinchPos, double _sigPenetFacPos, double _bilinEndAmpPos, double _unloadingStiffFacPos,
	int _cyclicRuleN, double _gapN, double _alphaPinchNeg, double _pinchYNeg,
	double _betaPinchNeg, double _epsPinchNeg, double _sigPenetFacNeg, double _bilinEndAmpNeg, double _unloadingStiffFacNeg,
	double sigInit) :
	UniaxialMaterial(tag, MAT_TAG_SmoothIMK),
	pd(_pd), pf(_pf),
	nd(_nd), nf(_nf),
	cS(_cS), FailEnergS(_gamaS* pf[0] * pd[0]),
	cUnloadE(_cUE), FailEnergUnloadE(_gamaUE* pf[0] * pd[0]),
	r0(_r0), r1(_r1), r2(_r2),
	cyclicRuleP(_cyclicRuleP), gapP(_gapP), alphaPinchPos(_alphaPinchPos), pinchYPos(_pinchYPos),
	betaPinchPos(_betaPinchPos), epsPinchPos(_epsPinchPos), sigPenetFacP(_sigPenetFacPos), bilinEndAmpP(_bilinEndAmpPos), unloadingStiffFacPos(_unloadingStiffFacPos),
	cyclicRuleN(_cyclicRuleN), gapN(_gapN), alphaPinchNeg(_alphaPinchNeg), pinchYNeg(_pinchYNeg),
	betaPinchNeg(_betaPinchNeg), epsPinchNeg(_epsPinchNeg), sigPenetFacN(_sigPenetFacNeg), bilinEndAmpN(_bilinEndAmpNeg), unloadingStiffFacNeg(_unloadingStiffFacNeg),
	sigini(sigInit)
{
	FyIndP = 0; // Index of the first positive force in the backbone
	FyIndN = 0; // Index of the first negative force in the backbone
	FrIndP = -1; //default behavior is no Fr when no descending backbone is given
	FrIndN = -1;
	// Find the index of the Fc and Fr values in pf
	auto it = std::max_element(pf.begin(), pf.end());
	FcIndP = std::distance(pf.begin(), it);
	if (it == pf.end() - 1)
	{
		//no descending branch
		FrIndP = -1;
	}
	else {
		it = std::min_element(pf.begin() + 1, pf.end());
		FrIndP = std::distance(pf.begin(), it);
		if (it == pf.end() - 1) {
			if (fabs((pf[FrIndP] - pf[FrIndP - 1]) / pf[FrIndP]) > 0.015)
			{
				//flat end line not detected; du should be computed by extrapolating the end line
				double du = pd[FrIndP] + (pd[FrIndP] - pd[FrIndP - 1]) / (pf[FrIndP] - pf[FrIndP - 1]) * (0 - pf[FrIndP]);
				pd[FrIndP] = du;
				pf[FrIndP] = 0;
			}
			else {
				FrIndP--;
			}
		}
		else {
			double min = *it;
			it++;
			for (; it < pf.end(); it++)
			{
				if (fabs((*it - min) / min) > 0.01)
				{
					opserr << "SmoothIMK:: Error the positive backbone increase by more than 1 percent after descending: tag = " << tag << endln;
					exit(-1);
				}
				*it = 0.99 * min;	//apply small increase in post-residual strength to prevent zero slope
			}
		}
	}
	// Find the index of the Fc and Fr values in nf
	it = std::min_element(nf.begin(), nf.end());
	FcIndN = std::distance(nf.begin(), it);
	if (it == nf.end() - 1)
	{
		//no descending branch
		FrIndN = -1;
	}
	else {
		it = max_element(nf.begin() + 1, nf.end());
		FrIndN = std::distance(nf.begin(), it);
		if (it == nf.end() - 1) {
			if (fabs((nf[FrIndN] - nf[FrIndN - 1]) / nf[FrIndN]) > 0.015)
			{
				//flat end line not detected; du should be computed by extrapolating the end line
				double du = nd[FrIndN] + (nd[FrIndN] - nd[FrIndN - 1]) / (nf[FrIndN] - nf[FrIndN - 1]) * (0 - nf[FrIndN]);
				nd[FrIndN] = du;
				nf[FrIndN] = 0;
			}
			else {
				FrIndN--;
			}
		}
		else {
			double max = *it;
			it++;
			for (; it < nf.end(); it++)
			{
				if (fabs((*it - max) / max) > 0.01)
				{
					opserr << "SmoothIMK:: Error (tag = " << tag << " ) the negative backbone increase by more than 1 percent after descending!" << endln;
					exit(-1);
				}
				*it = 0.99 * max;	//apply small increase in post-residual strength to prevent zero slope
			}
		}
	}
	revertToStart();
}

SmoothIMK::SmoothIMK(void) :
	UniaxialMaterial(0, MAT_TAG_SmoothIMK)
{
}

SmoothIMK::~SmoothIMK(void)
{
	// Does nothing
}

UniaxialMaterial*
SmoothIMK::getCopy(void)
{
	SmoothIMK* theCopy = new SmoothIMK(this->getTag(), pd, pf, nd, nf, FailEnergS / pf[0] / pd[0], cS, FailEnergUnloadE / pf[0] / pd[0],
		cUnloadE, r0, r1, r2,
		cyclicRuleP, gapP, alphaPinchPos, pinchYPos, betaPinchPos, epsPinchPos, sigPenetFacP, bilinEndAmpP, unloadingStiffFacPos,
		cyclicRuleN, gapN, alphaPinchNeg, pinchYNeg, betaPinchNeg, epsPinchNeg, sigPenetFacN, bilinEndAmpN, unloadingStiffFacNeg,
		sigini);
	theCopy->revertToStart();
	return theCopy;
}

double
SmoothIMK::getInitialTangent(void)
{
	return E0p;
}

double SmoothIMK::getStrain(void)
{
	return eps;
}

double SmoothIMK::getStress(void)
{
	return sig;
}

double SmoothIMK::getTangent(void)
{
	return e;
}

int
SmoothIMK::setParameter(const char** argv, int argc, Parameter& param)
{
	return -1;
}

int
SmoothIMK::updateParameter(int parameterID, Information& info)
{
	return -1;
}

int
SmoothIMK::commitState(void)
{
	epsminP = epsmin;
	epsmaxP = epsmax;
	epsLimitP = epsLimit;
	epss0P = epss0;
	sigs0P = sigs0;
	epssrP = epsr;
	sigsrP = sigr;
	slopeRatP = slopeRat;
	onEnvelopeP = onEnvelope;
	epsPlP = epsPl;
	updateDamage();
	isPosDirP = isPosDir;
	branchP = branch;
	eP = e;
	sigP = sig;
	epsP = eps;
	R0P = R0;
	initiatedP = initiated;
	return 0;
}

int
SmoothIMK::revertToLastCommit(void)
{
	epsmin = epsminP;
	epsmax = epsmaxP;
	epsLimit = epsLimitP;
	epss0 = epss0P;
	sigs0 = sigs0P;
	epsr = epssrP;
	sigr = sigsrP;
	branch = branchP;
	isPosDir = isPosDirP;
	slopeRat = slopeRatP;
	R0 = R0P;
	e = eP;
	sig = sigP;
	eps = epsP;
	onEnvelope = onEnvelopeP;
	epsPl = epsPlP;
	initiated = initiatedP;
	return 0;
}

int
SmoothIMK::revertToStart(void)
{
	EnergyP = 0;	//by SAJalali
	fpDmgd = pf;
	fnDmgd = nf;
	E0p = pf[0] / pd[0];
	E0n = nf[0] / nd[0];
	e = eP = E0p;
	onEnvelope = onEnvelopeP = true;
	branchP = branch = 0; //precap
	epsP = sigini / E0p;
	sigP = sigini;
	sig = 0.0;
	eps = 0.0;
	EunloadP = E0p * unloadingStiffFacPos;
	EunloadN = E0n * unloadingStiffFacNeg;
	isPosDirP = isPosDir = true;
	//epsmaxP = 0;
	//epsminP = 0;
	epsmaxP = pd[0];
	epsminP = nd[0];
	epsLimitP = pd[0];
	epss0P = 0.0;
	sigs0P = 0.0;
	epssrP = 0.0;
	sigsrP = 0.0;
	R0P = r0;
	epsPl = epsPlP = 0;
	ExcurEnergy = 0;
	slopeRat = slopeRatP = 0;
	initiated = initiatedP = false;
	// Copy pd to gpd and nd to gnd
	gpd = pd;
	gnd = nd;
	if (gapP > 0.001 * pd[0] || gapN < 0.001 * nd[0])
	{
		branchP = 10010; // gapped loading;
		for (size_t i = 0; i < pd.size(); ++i) {
			gpd[i] += gapP;
		}
		epsmax = gpd[0];
		for (size_t i = 0; i < nd.size(); ++i) {
			gnd[i] += gapN;
		}
		epsmin = gnd[0];
	}

	return 0;
}

int
SmoothIMK::sendSelf(int commitTag, Channel& theChannel)
{
	int numPData = pd.size();
	int numNData = nd.size();
	static Vector data(49 + (numPData + numNData) * 2);// 47 is the number of fixed data members
	int n = -1;
	data(n++) = numPData;	// size of the data vector
	data(n++) = numNData;	// size of the data vector
	data(n++) = this->getTag();	//0
	for (size_t i = 0; i < pd.size(); ++i) {
		data(n++) = pd[i];
		data(n++) = pf[i];
	}
	for (size_t i = 0; i < nd.size(); ++i) {
		data(n++) = nd[i];
		data(n++) = nf[i];
	}
	data(n++) = alphaPinchPos;
	data(n++) = pinchYPos;
	data(n++) = betaPinchPos;
	data(n++) = epsPinchPos;
	data(n++) = sigPenetFacP;
	data(n++) = alphaPinchNeg;
	data(n++) = pinchYNeg;
	data(n++) = betaPinchNeg;
	data(n++) = epsPinchNeg;
	data(n++) = sigPenetFacN;
	data(n++) = FailEnergS;
	data(n++) = cS;
	data(n++) = FailEnergUnloadE;
	data(n++) = cUnloadE;
	data(n++) = sigini;
	data(n++) = epsP;
	data(n++) = sigP;
	data(n++) = eP;
	data(n++) = EnergyP;
	data(n++) = epsminP;
	data(n++) = epsmaxP;
	data(n++) = epsLimitP;
	data(n++) = epss0P;
	data(n++) = sigs0P;
	data(n++) = epssrP;
	data(n++) = sigsrP;
	data(n++) = isPosDirP;
	data(n++) = branchP;
	data(n++) = FyIndP;
	data(n++) = FyIndN;
	data(n++) = FcIndP;
	data(n++) = FrIndP;
	data(n++) = FcIndN;
	data(n++) = FrIndN;
	data(n++) = epsPlP;
	data(n++) = ExcurEnergy;
	data(n++) = slopeRatP;
	data(n++) = onEnvelopeP;
	data(n++) = initiatedP;
	data(n++) = R0P;
	data(n++) = gapP;
	data(n++) = gapN;
	data(n++) = bilinEndAmpP;
	data(n++) = bilinEndAmpN;
	data(n++) = cyclicRuleP;
	data(n++) = cyclicRuleN;

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int
SmoothIMK::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	static Vector data(69);	//assuming 45 fixed data members plus 20 positive and negative backbone data

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "SmoothIMK::recvSelf() - failed to recvSelf\n";
		return -1;
	}
	int n = -1;
	int numPData = data(n++);	// size of the data vector
	int numNData = data(n++);	// size of the data vector
	if (numPData + numNData > 20) {
		opserr << "SmoothIMK::recvSelf() - Error: the number of positive and negative data exceeded the transfer limit\n";
		return -1;
	}
	this->setTag(data(n++));
	// Check the size of the data vector
	pd.resize(numPData);
	pf.resize(numPData);
	nd.resize(numNData);
	nf.resize(numNData);
	for (size_t i = 0; i < numPData; ++i) {
		pd.push_back(data(n++));
		pf.push_back(data(n++));
	}
	for (size_t i = 0; i < numNData; ++i) {
		nd.push_back(data(n++));
		nf.push_back(data(n++));
	}
	alphaPinchPos = data(n++);
	pinchYPos = data(n++);
	betaPinchPos = data(n++);
	epsPinchPos = data(n++);
	sigPenetFacP = data(n++);
	alphaPinchNeg = data(n++);
	pinchYNeg = data(n++);
	betaPinchNeg = data(n++);
	epsPinchNeg = data(n++);
	sigPenetFacN = data(n++);
	FailEnergS = data(n++);
	cS = data(n++);
	FailEnergUnloadE = data(n++);
	cUnloadE = data(n++);
	sigini = data(n++);
	epsP = data(n++);
	sigP = data(n++);
	eP = data(n++);
	EnergyP = data(n++);
	epsminP = data(n++);
	epsmaxP = data(n++);
	epsLimitP = data(n++);
	epss0P = data(n++);
	sigs0P = data(n++);
	epssrP = data(n++);
	sigsrP = data(n++);
	isPosDirP = data(n++);
	branchP = data(n++);
	FyIndP = data(n++);
	FyIndN = data(n++);
	FcIndP = data(n++);
	FrIndP = data(n++);
	FcIndN = data(n++);
	FrIndN = data(n++);
	epsPlP = data(n++);
	ExcurEnergy = data(n++);
	slopeRatP = data(n++);
	onEnvelopeP = data(n++);
	initiatedP = data(n++);
	R0P = data(n++);
	gapP = data(n++);
	gapN = data(n++);
	bilinEndAmpP = data(n++);
	bilinEndAmpN = data(n++);
	cyclicRuleP = data(n++);
	cyclicRuleN = data(n++);
	e = eP;
	sig = sigP;
	eps = epsP;
	//TODO: fill gpd and gnd based on the gap value and the initial loading direction
	return 0;
}

void
SmoothIMK::Print(OPS_Stream& s, int flag)
{
	//    s << "SmoothIMK:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
	const char* endStr = endln;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
	{
		endStr = "";
		s << "\t\t\t{";
	}
	s << "\"name \": \"" << this->getTag() << "\", " << endStr;
	s << "\"type\": \"SmoothIMK\", " << endStr;
	//s << "\"pd1 \": " << pd1 << ", " << endStr;
	//s << "\"pd2 \": " << pd2 << ", " << endStr;
	//s << "\"pd3 \": " << pd3 << ", " << endStr;
	//s << "\"nd1 \": " << nd1 << ", " << endStr;
	//s << "\"nd2 \": " << nd2 << ", " << endStr;
	//s << "\"nd3 \": " << nd3 << ", " << endStr;
	//s << "\"pf1 \": " << pf1 << ", " << endStr;
	//s << "\"pf2 \": " << pf2 << ", " << endStr;
	//s << "\"pf3 \": " << pf3 << ", " << endStr;
	//s << "\"nf[0] \": " << nf[0] << ", " << endStr;
	//s << "\"nf2 \": " << nf2 << ", " << endStr;
	//s << "\"nf3 \": " << nf3 << ", " << endStr;
	s << "\"sigini \": " << sigini << endStr;
	if (flag == OPS_PRINT_PRINTMODEL_JSON)
		s << "}";

}

Response* SmoothIMK::setResponse(const char** argv, int argc, OPS_Stream* theOutput)
{
	Response* theResponse = UniaxialMaterial::setResponse(argv, argc, theOutput);
	if (theResponse != 0)
		return theResponse;
	if (strcmp(*argv, "branch") == 0)
	{
		theResponse = new MaterialResponse(this, 11, 0);
		return theResponse;
	}
	return 0;
}

int SmoothIMK::getResponse(int responseID, Information& matInformation)
{
	int res = 0;
	res = UniaxialMaterial::getResponse(responseID, matInformation);
	if (res == 0)
		return 0;
	switch (responseID)
	{
	case 11:
		matInformation.setInt(branchP);
		return 0;
	default:
		return -1;
	}
}

void SmoothIMK::updateDamage()
{
	if (sigP * sig < 0)
	{
		if (sigP > 0)
		{
			double zeroSigEps = epsP - sigP / E0p;
			double dE = 0.5 * sigP * (zeroSigEps - epsP);
			EnergyP += dE;
			if (EnergyP < 0) EnergyP = 0.;
			ExcurEnergy += dE;
			if (ExcurEnergy < 0) ExcurEnergy = 0.;
			if (branch == pd.size() - 1)
			{
				// we are on residual branch
				FyIndP = FrIndP;
				FcIndP = FrIndP;
			}
			else if (branch < 1000 && branch >= pd.size())
			{
				FyIndP = -1;
				FcIndP = -1;
			}
			else
			{
				double beta = 0;
				if (EnergyP > FailEnergS)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergS - EnergyP), cS);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Strength loss\n" << endln;
					beta = 0.9999;
				}
				double Fr = FrIndP != -1 ? pf[FrIndP] : 0.00001*pf[0];
				double Fr_ = 0;
				for (size_t i = 0; i < FrIndP; ++i) {
					Fr_ = pf[i] / pf[0] * Fr;
					fpDmgd[i] = (1. - beta) * (fpDmgd[i] - Fr_) + Fr_;
				}
				//FrP = FcP - pf2 + pf3;
				//if (FrP < 0)
				//	FrP = 0;

				if (EnergyP > FailEnergUnloadE)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergUnloadE - EnergyP), cUnloadE);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Unloading Stiffness loss\n" << endln;
					beta = 0.9999;
				}
				EunloadP *= (1. - beta);
			}
		}
		else //sigP < 0
		{
			double zeroSigEps = epsP - sigP / E0n;
			double dE = 0.5 * sigP * (zeroSigEps - epsP);
			EnergyP += dE;
			if (EnergyP < 0) EnergyP = 0.;
			ExcurEnergy += dE;
			if (ExcurEnergy < 0) ExcurEnergy = 0.;
			if (branch == pd.size() - 1)
			{
				FyIndN = FrIndN;
				FcIndN = FrIndN;
			}
			else if (branch < 1000 && branch >= pd.size())
			{
				FyIndN = -1;
				FcIndN = -1;
			}
			else {
				double beta = 0;
				if (EnergyP > FailEnergS)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergS - EnergyP), cS);
				if (beta > 0.999 || beta < 0)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Strength loss\n" << endln;
					beta = 0.999;
				}
				double Fr = FrIndN != -1 ? nf[FrIndN] : 0.00001 * nf[0];
				double Fr_ = 0;
				for (size_t i = 0; i < FrIndN; ++i) {
					Fr_ = nf[i] / nf[0] * Fr;
					fnDmgd[i] = (1. - beta) * (fnDmgd[i] - Fr_) + Fr_;
				}
				//FrN = FcN - nf2 + nf3;
				//if (FrN > 0)
				//	FrN = 0;
				if (EnergyP > FailEnergUnloadE)
					beta = 1;
				else
					beta = pow(ExcurEnergy / (FailEnergUnloadE - EnergyP), cUnloadE);
				if (beta > 0.9999)
				{
					opserr << "\nSmoothIMK:" << this->getTag() << " WARNING! Complete Unloading Stiffness loss\n" << endln;
					beta = 0.9999;
				}
				EunloadN *= (1. - beta);
			}
		}
		ExcurEnergy = 0.0;
	}
	else
	{
		double dE = 0.5 * (sig + sigP) * (eps - epsP);
		ExcurEnergy += dE;
		EnergyP += dE;
	}
}

int SmoothIMK::setTrialStrain(double trialStrain, double strainRate)
{
	revertToLastCommit();
	eps = trialStrain + (sigini > 0 ? sigini / E0p : sigini / E0n);
	const double deps = eps - epsP;
	if (fabs(deps) < 1.e-20)
		return 0;
	if (!initiated) // handle start in negative direction
	{
		isPosDir = (deps > 0);
		initiated = true;
		updateAsymptote();
	}
	else if (isPosDir && deps < 0.0)
	{
		isPosDir = false;
		if (sigP >= 0)
		{
			epsPl = epsP - sigP / EunloadP;
		}
		if (epsP > epsmax)
			epsmax = epsP;
		onEnvelope = false;
		changeBranch(sigP < 0);
		updateAsymptote();
	}
	else if (deps > 0.0 && !isPosDir) {
		isPosDir = true;
		if (sigP <= 0)
		{
			epsPl = epsP - sigP / EunloadN;
		}
		if (epsP < epsmin)
			epsmin = epsP;
		onEnvelope = false;
		changeBranch(sigP > 0);
		updateAsymptote();
	}
	else if ((isPosDir && eps > epsLimit) || (!isPosDir && eps < epsLimit))
	{
		onEnvelope = true;
		changeBranch(false);
		updateAsymptote();
	}
	double epsrat = (eps - epsr) / (epss0 - epsr);
	double dum1 = 1.0 + pow(fabs(epsrat), R0);
	double dum2 = pow(dum1, (1 / R0));

	sig = slopeRat * epsrat + (1.0 - slopeRat) * epsrat / dum2;
	sig = sig * (sigs0 - sigr) + sigr;

	e = slopeRat + (1.0 - slopeRat) / (dum1 * dum2);
	e *= (sigs0 - sigr) / (epss0 - epsr);
	return 0;
}

void SmoothIMK::changeBranch(bool isReturning)
{
	bool shouldPinch = (
		(isPosDir && fabs(epsPl - epsmax) / pd[0] >= bilinEndAmpP) ||
		(!isPosDir && fabs(epsPl - epsmin) / nd[0] <= -bilinEndAmpN)
		);
	if (onEnvelope)
	{
		branch = nextBranch(branch, shouldPinch);
		return;
	}
	//unloading:
	int& cyclicRule = isPosDir ? cyclicRuleP : cyclicRuleN;
	double& gap = isPosDir ? gapP : gapN;
	std::vector<double>& d = isPosDir ? pd : nd;
	bool hasGap = isPosDir ? gap > 0.001 * d[0] : gap < 0.001 * d[0];
	if (hasGap && !isReturning)
	{
		if (cyclicRule == 3 && shouldPinch)
			branch = 1003;
		else
			branch = 10012;
		return;
	}
	branch = cyclicRule == 1 || !shouldPinch ? 0 : cyclicRule == 2 ? 1006 : 1003;
}

int SmoothIMK::nextBranch(int branch, bool shouldPinch)
{
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const std::vector<double>& d = isPosDir ? pd : nd;
	int cyclicRule = isPosDir ? cyclicRuleP : cyclicRuleN;
	if (branch == 10010 || branch == 10011)
		return 0;
	else if (branch == 10012)
	{
		if (cyclicRule == 1 || !shouldPinch)
			return 10011;
		return 1002;
	}
	else if (branch == 1002 || branch == 1004 || branch == 1006)
		return 1005;
	else if (branch == 1003)
		return 1004;
	else if (branch == 1005)
	{
		double tmp;
		int targBranch;
		getEnvelope(epsPeak, false, tmp, targBranch, tmp, tmp);
		return targBranch;
	}
	else if (branch == 1010 || (branch < 1000 && branch > d.size() - 1))
		return 1010;//failed
	else
		return branch + 1;

}

void SmoothIMK::updateAsymptote()
{
	//updates: epsr, sigr, epss0, sigs0, epsLimit, slopeRat, R0
	const std::vector<double>& gd = isPosDir ? gpd : gnd;
	const std::vector<double>& d = isPosDir ? pd : nd;
	const std::vector<double>& f = isPosDir ? fpDmgd : fnDmgd;
	const int& FyInd = isPosDir ? FyIndP : FyIndN;
	const int& FcInd = isPosDir ? FcIndP : FcIndN;
	const int& FrInd = isPosDir ? FrIndP : FrIndN;
	const double& E1 = isPosDir ? E0p : E0n;
	const double Esh = FyInd == -1 ? 0.001*E1 : (f[FyInd + 1] - f[FyInd]) / (gd[FyInd + 1] - gd[FyInd]);
	const double& Fy = FyInd == -1 ? 0.001*f[0] : f[FyInd];
	const double& dy = gd[0];
	const double& dc = gd[FcInd];
	const double& E2 = isPosDir ? (sigP < -0.001 * Fy ? EunloadN : E0p) : (sigP > 0.001 * Fy ? EunloadP : E0n);
	const double& Eunload = isPosDir ? EunloadN : EunloadP;
	const double& gapD = isPosDir ? gapP : gapN;
	const double  gap = (gapP - gapN)/2; // sum abs values by 2
	const double& epsPeak = isPosDir ? epsmax : epsmin;
	const double& sigPenetFac = isPosDir ? sigPenetFacP : sigPenetFacN;
	const double& alphaPinch = isPosDir ? alphaPinchPos : alphaPinchNeg;
	const double& pinchY = isPosDir ? pinchYPos : pinchYNeg;
	const double& betaPinch = isPosDir ? betaPinchPos : betaPinchNeg;
	const double& epsPinch = isPosDir ? epsPinchPos : epsPinchNeg;
	bool shouldPinch = (
		(isPosDir && fabs(epsPl - epsmax) / pd[0] >= bilinEndAmpP) ||
		(!isPosDir && fabs(epsPl - epsmin) / nd[0] <= -bilinEndAmpN)
		);
	epsr = epsP;
	sigr = sigP;
	double k2 = 0;
	double Dy = 0;
	if (branch == 10010 || branch == 10011) // epsPl+gapD to yield
	{
		//gapped loading
		epss0 = epsPl + (branch == 10010 ? gapD : gap * 2);
		sigs0 = (branch == 10010 ? 0.001 : 0.002) * Fy;
		Dy = (Fy - Esh * dy - sigs0 + E2 * epss0) / (E2 - Esh);
		epsLimit = (epss0 + Dy) / 2;
		k2 = E1;
	}
	else if (branch == 10012) // epsPl to epsPl+gapD
	{
		//gapped unloading
		sigs0 = 0.00001 * Fy;
		if (fabs(epsPl - epsr) < 0.01 * fabs(d[0]))
		{ // unloading from gap should be skipped
			branch = nextBranch(10012, shouldPinch);
			return updateAsymptote();
		}
		epss0 = epsPl;
		epsLimit = epss0 + gap;
		k2 = 0.001 * Fy / gap;
	}
	else if (branch == 0) // yield to 1
	{
		epss0 = (Fy - Esh * dy - sigr + E2 * epsr) / (E2 - Esh);
		if ((!isPosDir && (epss0 > epsr || epss0 < dc)) || (isPosDir && (epss0 < epsr || epss0 > dc)))
		{
			//we should skip this branch
			branch++;
			return updateAsymptote();
		}
		sigs0 = (epss0 - (epsP - sigP / E2)) * E1;
		epsLimit = epss0 + (gd[1] - gd[0]) / 2;
		k2 = Esh;
	}
	else if (branch == 1002 || branch == 1006) // epsPl+gapD to peak
	{
		if (branch == 1002)
		{
			epss0 = epsPl + 2 * gap;
			sigs0 = 0.002 * Fy;
		}
		else { //1006
			epss0 = epsPl;
			sigs0 = (sigPenetFac > 0.00001 ? sigPenetFac : 0.00001) * Fy;
		}
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		int targBranch;
		double tmp, nextLimit, nextEps = epsPeak;
		getEnvelope(epsPeak, false, sigmax, targBranch, tmp, nextLimit);
		if (fabs((epsPeak - nextLimit) / d[0]) < 0.1	)
		{
			targBranch = nextBranch(1005, shouldPinch);
			if (targBranch < gd.size())
			{
				nextEps = gd[targBranch];
				sigmax = f[targBranch];
			}
			else {
				nextEps = 1000;
				sigmax = 0;
			}
		}
		k2 = (sigmax - sigs0) / (nextEps - epss0);
	}
	else if (branch == 1003) // epsPl to epsPl+gapD+alphaPinch
	{
		double sigmax;
		int targBranch; // not used
		double tmp;
		getEnvelope(epsPeak, false, sigmax, targBranch, tmp, tmp);
		double penetFac = (sigPenetFac > 0.00001 ? sigPenetFac : 0.00001);
		sigs0 = penetFac *Fy;
		epss0 = epsPl + penetFac * d[0];
		double x = epsPeak - epss0 - 2*gap;
		double dx = pinchDx(x, alphaPinch, epsPinch, betaPinch, d[0]);
		double pinchEps = epss0 + dx + 2*gap;
		epsLimit = (epss0 + pinchEps) / 2;
		double y = pinchY * sigmax;
		k2 = (y - sigs0) / (pinchEps - epss0);
	}
	else if (branch == 1004) // epsPl+gapD+alphaPinch to peak
	{
		double penetFac = (sigPenetFac > 0.00001 ? sigPenetFac : 0.00001);
		double x = epsPeak - epsPl - penetFac * d[0] - 2*gap;
		double dx = pinchDx(x, alphaPinch, epsPinch, betaPinch, d[0]);
		epss0 = epsPl + penetFac * d[0] + dx + 2*gap;
		epsLimit = (epss0 + epsPeak) / 2;
		double sigmax;
		int targBranch;
		double tmp, nextLimit, nextEps = epsPeak;
		getEnvelope(epsPeak, false, sigmax, targBranch, tmp, nextLimit);
		sigs0 = pinchY * sigmax;
		if (fabs((epsPeak - nextLimit) / d[0]) < 0.1)
		{
			targBranch = nextBranch(1005, shouldPinch);
			if (targBranch < gd.size())
			{
				nextEps = gd[targBranch];
				sigmax = f[targBranch];
			}
			else {
				nextEps = 1000;
				sigmax = 0;
			}
		}
		k2 = (sigmax - sigs0) / (nextEps - epss0);
	}
	else if (branch == 1005)// peak-passing
	{
		epss0 = epsPeak;
		int targBranch;
		getEnvelope(epss0, true, sigs0, targBranch, k2, epsLimit);
		if (fabs((epss0 - epsLimit) / d[0]) < 0.1)
		{
			branch = targBranch;
			return updateAsymptote();
		}
	}
	else if (branch == 1010) //failed branch
	{

		int b = d.size();
		const double& du = d[b - 1];
		const double& dr = d[b - 2];
		epss0 = (2 * du - dr);
		sigs0 = 0.00001*Fy;
		epsLimit = d[0] * 1000;
		k2 = 1e-5*E1;
	}
	else if (branch >= 1 && branch < 1000)
	{
		if (branch < d.size() - 1)
		{
			epss0 = gd[branch];
			sigs0 = f[branch];
			double d2, f2;
			d2 = gd[branch + 1];
			f2 = f[branch + 1];
			k2 = (f2 - sigs0) / (d2 - epss0);
			epsLimit = (epss0 + d2) / 2;
		}
		else if (branch < d.size())
		{
			//failing branch
			epss0 = gd[branch];
			sigs0 = f[branch];
			const double& du = gd[branch];
			const double& dr = gd[branch - 1];
			epsLimit = (3 * du - dr) / 2;
			k2 = (f[branch - 1] - f[branch - 2]) / (gd[branch - 1] - gd[branch - 2]);
			Dy = du - f[branch] / k2;
			if ((isPosDir && Dy < epsLimit) || (!isPosDir && Dy > epsLimit))
				epsLimit = Dy;
		}
	}
	else {
		opserr << "ERROR in SmoothIMK::updateAsymptote: unrecognized branch: " << branch << endln;
		exit(-1);
	}
	double k1 = (sigs0 - sigr) / (epss0 - epsr);
	slopeRat = k2 / k1;
	computeR0(k1, k2, E1, dy);
}

void SmoothIMK::getEnvelope(double eps, bool onPlusSide, double& targStress, int& targBranch, double& k, double& limitEps)
{
	const std::vector<double>& dVec = eps > 0 ? gpd : gnd;
	const std::vector<double>& fVec = eps > 0 ? fpDmgd : fnDmgd;

	// Interpolate stress based on fabs(eps), dVec, fVec
	size_t n = dVec.size();
	if (fabs(eps) >= fabs(dVec.back())) {
		// Extrapolate above maximum
		k = (fVec.back() - fVec[n - 2]) / (dVec.back() - dVec[n - 2]);
		targStress = fVec.back() + (eps - dVec.back()) * k;
		limitEps = dVec.back() * 100;
		targBranch = dVec.size() - 1;
	}
	else {
		// Interpolate between points
		for (size_t i = 0; i < n; ++i) {
			if (fabs(eps) < fabs(dVec[i]) || fabs((eps - dVec[i]) / dVec[0]) < 0.01) {
				targBranch = i;
				double k1 = (fVec[i] - fVec[i - 1]) / (dVec[i] - dVec[i - 1]);
				k = k1;
				if (fabs((eps - dVec[i]) / dVec[0]) < 0.01) {
					//we are on a vertex and need to decide about which side
					if (i < n)
						limitEps = (dVec[i] + dVec[i + 1]) / 2;
					else
						limitEps = 2 * dVec[i];
					if (onPlusSide && i < n - 1)
						k = (fVec[i + 1] - fVec[i]) / (dVec[i + 1] - dVec[i]);
					targStress = fVec[i];
				}
				else
				{
					double d1 = epss0;
					if (i > 0)
						d1 = dVec[i - 1];
					double tmp = (d1 + dVec[i]) / 2;
					if (fabs(eps) < fabs(tmp))
						limitEps = tmp;
					else
						limitEps = (dVec[i] + eps) / 2;
					targStress = fVec[i - 1] + (eps - dVec[i - 1]) * k1;
				}
				break;
			}
		}
	}
}

void SmoothIMK::computeR0(double k1, double k2, double E1, double dy)
{
	R0 = r0;
	if (r1 != 0.0)
	{
		double xi_1 = E1 / fabs(k1 - k2);
		double xi_21 = fabs(epss0 - epsr) / dy;
		double xi_22 = fabs(epsLimit - epss0) / dy;
		double xi_2 = (2 * std::min(xi_21, xi_22));
		R0 += r1 * xi_1 + r2 * xi_2;
		R0 = std::min(30.0, R0);
		R0 = std::max(1.0, R0);
	}
}
