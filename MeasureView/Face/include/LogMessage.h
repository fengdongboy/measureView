#pragma once
#ifndef __LOGMESSAGE_H__
#define __LOGMESSAGE_H__

#include <string>

#define log1(a) LogMessage::getSingletonPtr()->logMessage(a)
#define log2(a, b) LogMessage::getSingletonPtr()->logMessage(a, b)
#define logv LogMessage::getSingletonPtr()->logMessagev

class LogMessage
{
public:
	LogMessage(void);
	~LogMessage(void);

	static LogMessage* getSingletonPtr( void );

	/// һ����Ϣ��
	void logMessage( const std::string& str );

	/// ������Ϣ��
	void logMessage( const std::string& str1, const std::string& str2);

	void logMessagev( const char* buf, ...);

	/// ����ļ�����
	void setLogName( const std::string& name )
	{
		mName = name;
	}

protected:
	void clear();

private:
	static LogMessage* mLogMessage;
	std::string mName;
};

#endif