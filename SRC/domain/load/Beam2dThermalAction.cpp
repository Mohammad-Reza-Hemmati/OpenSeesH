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

// $Revision: 1.1 $
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dThermalAction.cpp,v $


//Modified by Jian Zhang, [Univeristy of Edinburgh]
//Modified by Panagiotis Kotsovinos, [Univeristy of Edinburgh]
//Modified by Liming Jiang [http://openseesforfire.github.io]

// Description: This file contains the class implementation for Beam2dThermalAction.
// Beam2dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Beam2dThermalAction.h>
#include <Vector.h>
#include <Element.h>
#include <Domain.h>
Vector Beam2dThermalAction::data(0);

Beam2dThermalAction::Beam2dThermalAction(int tag,
	double t1, double locY1, double t2, double locY2,
	double t3, double locY3, double t4, double locY4,
	double t5, double locY5, double t6, double locY6,
	double t7, double locY7, double t8, double locY8,
	double t9, double locY9,
	int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag),
	ThermalActionType(LOAD_TAG_Beam2dThermalAction), theSeries(0)
{
	Temp.resize(9);
	Temp[0] = t1; Temp[1] = t2; Temp[2] = t3; Temp[3] = t4; Temp[4] = t5;
	Temp[5] = t6; Temp[6] = t7; Temp[7] = t8; Temp[8] = t9;
	Loc[0] = locY1; Loc[1] = locY2; Loc[2] = locY3; Loc[3] = locY4; Loc[4] = locY5;
	Loc[5] = locY6; Loc[6] = locY7; Loc[7] = locY8; Loc[8] = locY9;

	Factors.Zero();
	indicator = 1; //without path timeseries defined;
}


Beam2dThermalAction::Beam2dThermalAction(int tag,
	double locY1, double locY2,
	TimeSeries* theSeries, int theElementTag
)
	:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), theSeries(theSeries),
	ThermalActionType(LOAD_TAG_Beam2dThermalAction)
{
	Loc.resize(9);
	Loc[0] = locY1;
	Loc[8] = locY2;

	for (int i = 1; i < 8; i++) {
		Loc[i] = Loc[0] - i * (Loc[0] - Loc[8]) / 8;
	}


	for (int i = 0; i < 9; i++) {
		Temp[i] = 0;
		TempApp[i] = 0;
	}
	Factors.Zero();
	indicator = 2;// Independent timeseries were created;

}

Beam2dThermalAction::Beam2dThermalAction(int tag,
	const Vector& locs,
	TimeSeries* theSeries, int theElementTag
)
	:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), theSeries(theSeries),
	ThermalActionType(LOAD_TAG_Beam2dThermalAction)
{

	for (int i = 0; i < locs.Size(); i++) {
		Loc[i] = locs(i);
	}

	Temp.resize(locs.Size());
	TempApp.resize(locs.Size());
	for (int i = 0; i < locs.Size(); i++) {
		Temp[i] = 0;
		TempApp[i] = 0;
	}
	Factors.Zero();
	indicator = 2;// Independent timeseries were created;

}

Beam2dThermalAction::Beam2dThermalAction(int tag,
	int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), ThermalActionType(LOAD_TAG_NodalThermalAction), theSeries(0)
{
	for (int i = 0; i < Temp.Size(); i++) {
		Temp[i] = 0;
		TempApp[i] = 0;
		Loc[i] = 0;
	}
	Factors.Zero();
	indicator = 3;// USing Nodal Thermal Action;
}

Beam2dThermalAction::Beam2dThermalAction()
	:ElementalLoad(LOAD_TAG_Beam2dThermalAction), ThermalActionType(LOAD_TAG_NodalThermalAction), theSeries(0)
{
	indicator = 3;// USing Nodal Thermal Action;
}

Beam2dThermalAction::~Beam2dThermalAction()
{
	indicator = 0;
	//if(theSeries!=0)
	//delete theSeries;

	theSeries = 0;
}

const Vector&
Beam2dThermalAction::getData(int& type, double loadFactor)
{
	type = ThermalActionType;
	if (data.Size() != 2 * TempApp.Size())
		data.resize(2 * TempApp.Size());
	for (int i = 0; i < TempApp.Size(); i++) {
		data(2 * i) = TempApp[i];
		data(2 * i + 1) = Loc[i];
	}
	Factors.Zero();
	return data;
}

void
Beam2dThermalAction::applyLoad(const Vector& factors)
{
	for (int i = 0; i < Temp.Size(); i++) {
		TempApp[i] = Temp[i] * factors(i);
	}
	if (theElement != 0)
		theElement->addLoad(this, factors(0));
}

//void Beam2dThermalAction::setDomain(Domain* pDomain)
//{
//	theElement = pDomain->getElement(eleTag);
//	if (theElement == 0)
//	{
//		opserr << "Beam2dThermalAction::setDomain() - failed to retrieve element pointer" << endln;
//	}
//}
//
void
Beam2dThermalAction::applyLoad(double loadfactor)
{
	// first determine the load factor
	if (indicator == 2 && theSeries != 0) {
		// theSeries == 0 occurs after recieveSelf where constant tempApp will be used
		Factors = ((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		for (int i = 0; i < TempApp.Size(); i++) {
			//PathTimeSeriesThermal returns absolute temperature;
			TempApp[i] = Factors(i);
		}
	}
	else if (indicator == 1) {
		for (int i = 0; i < Temp.Size(); i++) {
			TempApp[i] = Temp[i] * loadfactor;
		}
	}

	if (theElement != 0)
		theElement->addLoad(this, loadfactor);
}

int
Beam2dThermalAction::sendSelf(int commitTag, Channel& theChannel)
{
	// Currently implemented to only support saving after thermomechanical analysis is complete
	// Therefore, the next analysis is assumed to constantly use the temp value reached at the end of previou phase.
	// So, no timeSeries is saved and restored.
	// Also, only indicator == 2 (timeSeries-dependant temperature) is supported.

	ID idData(5);
	idData(0) = this->getTag();
	idData(1) = eleTag;
	idData(2) = indicator;
	idData(3) = ThermalActionType;
	idData(4) = TempApp.Size();
	int res = theChannel.sendID(this->getDbTag(), commitTag, idData);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::sendSelf() - failed to send Id data" << endln;
		return res;
	}
	opserr << "sent tempApp:\n" << TempApp << endln;
	res = theChannel.sendVector(this->getDbTag(), commitTag, TempApp);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::sendSelf() - failed to send TempApp" << endln;
		return res;
	}
	res = theChannel.sendVector(this->getDbTag(), commitTag, Loc);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::sendSelf() - failed to send Loc" << endln;
		return res;
	}
	return 0;
}

int
Beam2dThermalAction::recvSelf(int commitTag, Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	ID idData(5);
	int res = theChannel.recvID(this->getDbTag(), commitTag, idData);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::recvSelf() - failed to recieve data" << endln;
		return res;
	}
	this->setTag(idData(0));
	eleTag = idData(1);
	indicator = idData(2);
	if (indicator != 2)
	{
		opserr << "Beam2dThermalAction::recvSelf(): WARNING, Temp data are not sent and recieved and only indicator == 1 is going to work" << endln;
	}
	ThermalActionType = idData(3);
	theElement = 0;
	TempApp.resize(idData(4));
	Loc.resize(idData(4));
	res = theChannel.recvVector(this->getDbTag(), commitTag, TempApp);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::sendSelf() - failed to recieve TempApp" << endln;
		return res;
	}
	opserr << "rcvd tempApp:\n" << TempApp << endln;

	res = theChannel.recvVector(this->getDbTag(), commitTag, Loc);
	if (res < 0)
	{
		opserr << "Beam2dThermalAction::sendSelf() - failed to recieve Loc" << endln;
		return res;
	}
	theSeries = 0;
	Temp.resize(0);
	return 0;
}

// do it later
void
Beam2dThermalAction::Print(OPS_Stream& s, int flag)
{
	s << "Beam2dThermalAction - reference load : " << Temp[0] << " change  temp of bot\n";
	s << Temp[8] << " change  temp at top\n";
	s << "  element acted on: " << eleTag << endln;
}

