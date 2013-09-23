/*
 * readscene.h
 *
 *  Created on: Mar 6, 2010
 *      Author: vik
 */

#ifndef READSCENE_H_
#define READSCENE_H_

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include "camera.h"
#include "globals.h"
#include "objects.h"

float getTokenAsFloat (std::string inString, int whichToken);

void parseSceneFile (char *filnam);

#endif /* READSCENE_H_ */
