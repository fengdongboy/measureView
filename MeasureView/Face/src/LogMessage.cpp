#include "LogMessage.h"
#include "stdafx.h"

LogMessage* LogMessage::mLogMessage = NULL;

LogMessage::LogMessage(void)
	:mName("faceLog.log")
{
	//assert(mLogMessage == NULL);
	clear();
	mLogMessage = this;
}


LogMessage::~LogMessage(void)
{
}

LogMessage* LogMessage::getSingletonPtr(void)
{
	static LogMessage* l = new LogMessage;
	return l;
}

void LogMessage::clear()
{
	FILE* f = fopen(mName.c_str(), "w");
	if (f)
		fclose(f);
}

void LogMessage::logMessagev( const char* buf, ... )
{
	FILE* file = fopen(mName.c_str(), "a");
	if (!file)
	{
		printf("face create log is wrong\n");
		return;
	}
	va_list argptr;
	va_start( argptr, buf );	
	vfprintf(file, buf, argptr);
	va_end( argptr );
	fprintf(file, "\n");
	fclose(file);
}

/// һ����Ϣ��
void LogMessage::logMessage( const std::string& str )
{
	logMessagev(str.c_str());
}

/// ������Ϣ��
void LogMessage::logMessage( const std::string& str1, const std::string& str2)
{
	logMessagev(str1.c_str(), str2.c_str());
}