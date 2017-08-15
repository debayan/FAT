#pragma once

#include "Interaction.h"

class  InterSet
{
public:
	 InterSet();
	~InterSet();

	void load(QString filename);
	void save(QString filename);

	QVector<Interaction*>  interactions;
};

