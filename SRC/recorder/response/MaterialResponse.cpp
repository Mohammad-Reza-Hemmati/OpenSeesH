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

// $Revision: 1.4 $
// $Date: 2007-03-02 00:12:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/recorder/response/MaterialResponse.cpp,v $

// Written: MHS 
// Created: Oct 2000
//
// Description: This file contains the MaterialResponse class implementation

#include <MaterialResponse.h>
#include <Material.h>
#include <MaterialState.h>

Information&
MaterialResponse::getInformation(int respNum)
{
	if (respNum == 0) {
		return myInfo;
	}
	return infoList[respNum];
}
const Vector& MaterialResponse::getData(void)
{
	if (numMaterialas == 0)
		return myInfo.getData();

	int sizeAll = 0;
	if (respVec == 0) {
		for (int i = 0; i < numMaterialas; i++) {
			const Vector& matData = infoList[i].getData();
			sizeAll += matData.Size();
		}
		respVec = new Vector(sizeAll);
	}
	sizeAll = 0;
	for (int i = 0; i < numMaterialas; i++) {
		const Vector& matData = infoList[i].getData();
		int sizeMat = matData.Size();
		for (int j = 0; j < sizeMat; j++)
			(*respVec)(sizeAll + j) = matData(j);
		sizeAll += sizeMat;
	}
	return *respVec;
}
MaterialResponse::MaterialResponse(Material* mat, int id) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
	infoList[0] = Information();
}

MaterialResponse::MaterialResponse(Material* mat, int id, int val) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
	infoList[0] = Information(val);
}

MaterialResponse::MaterialResponse(Material* mat, int id, double val) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
	infoList[0] = Information(val);
}

MaterialResponse::MaterialResponse(Material* mat, int id, const ID& val) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
	infoList[0] = Information(val);
}

MaterialResponse::MaterialResponse(Material* mat, int id, const Vector& val) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
  new (&infoList[0]) Information(val);
}

MaterialResponse::MaterialResponse(Material* mat, int id, const Matrix& val) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(1), respVec(0)
{
	theMaterials = new Material * [1];
	theMaterials[0] = mat;
	infoList = new Information[1];
	new (&infoList[0]) Information(val);
}
MaterialResponse::MaterialResponse(Material** mats, int id, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		new (&infoList[i]) Information();
	}
}

MaterialResponse::MaterialResponse(Material** mats, int id, int val, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		new (&infoList[i]) Information(val);
	}
}

MaterialResponse::MaterialResponse(Material** mats, int id, double val, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		new (&infoList[i]) Information(val);
	}
}

MaterialResponse::MaterialResponse(Material** mats, int id, const ID& val, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		new (&infoList[i]) Information(val);
	}
}

MaterialResponse::MaterialResponse(Material** mats, int id, const Vector& val, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		new (&infoList[i]) Information(val);
	}
}

MaterialResponse::MaterialResponse(Material** mats, int id, const Matrix& val, int numMats) :
	Response(), theMaterialState(0), responseID(id), numMaterialas(numMats), respVec(0)
{
	theMaterials = new Material * [numMats];
	infoList = new Information[numMats];
	for (int i = 0; i < numMats; i++)
	{
		theMaterials[i] = mats[i];
		infoList[i] = Information(val);
	}
}

MaterialResponse::MaterialResponse(MaterialState* mat, int id) :
	Response(), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::MaterialResponse(MaterialState* mat, int id, int val) :
	Response(val), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::MaterialResponse(MaterialState* mat, int id, double val) :
	Response(val), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::MaterialResponse(MaterialState* mat, int id, const ID& val) :
	Response(val), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::MaterialResponse(MaterialState* mat, int id, const Vector& val) :
	Response(val), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::MaterialResponse(MaterialState* mat, int id, const Matrix& val) :
	Response(val), theMaterials(0), theMaterialState(mat), responseID(id), numMaterialas(0), respVec(0)
{

}

MaterialResponse::~MaterialResponse()
{
	for (int i = 0; i < numMaterialas; i++)
		theMaterials[i]->resetResponse(responseID, &myInfo);
	delete[] theMaterials;
	delete[] infoList;
	if (respVec != 0)
		delete respVec;
}

int
MaterialResponse::getResponse(void)
{
	if (theMaterials != 0)
	{
		int res = 0;
		for (int i = 0; i < numMaterialas; i++)
			res += theMaterials[i]->getResponse(responseID, infoList[i]);
		return res;
	}
	else if (theMaterialState != 0)
		return theMaterialState->getResponse(responseID, myInfo);
	else
		return 0;
}

int
MaterialResponse::getResponseSensitivity(int gradNumber)
{
	if (theMaterials != 0)
	{
		int res = 0;
		for (int i = 0; i < numMaterialas; i++)
			res += theMaterials[i]->getResponseSensitivity(responseID, gradNumber, infoList[i]);
		return res;
	}
}
