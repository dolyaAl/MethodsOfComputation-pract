#include "FunctionReader.h"
#include <algorithm>

void updateFunStr(std::string& str)
{
	str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
}
