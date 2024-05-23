#ifndef objmgr_h
#define objmgr_h

#include "objtype.h"
#include "obj.h"

#include "seqinfo.h"
#include "alignresult.h"
#include "pathinfo.h"

class Obj;

#define T(x)	class x;
#include "objtypes.h"

const char *ObjTypeToStr(ObjType Type);

class ObjMgr
	{
	friend class Obj;

public:
	Obj *m_Free[OTCount];
	Obj *m_Busy[OTCount];
#if DEBUG
	bool m_Validate;
#endif

private:
	ObjMgr();

public:
	static vector<ObjMgr *> m_OMs;
	static vector<ObjMgr *> m_FreeOMs;

public:
	static ObjMgr *CreateObjMgr();
	static void FreeObjMgr(ObjMgr *OM);

public:
	void Up(Obj *pObj);
	void Down(Obj *pObj);
	uint GetBusyCount(uint Type) const;
	uint GetFreeCount(uint Type) const;
	void LogStats() const;
#if TRACK_OBJS
	Obj *GetObj(ObjType Type, const char *FileName, uint LineNr);
#else
	Obj *GetObj(ObjType Type);
#endif

#if TRACK_OBJS
#define T(x)	\
	x *__Get##x(const char *FileName, uint LineNr) \
		{ return (x *) GetObj(OT_##x, FileName, LineNr); }
#define GetSeqInfo()		__GetSeqInfo(__FILE__, __LINE__)
#define GetPathInfo()		__GetPathInfo(__FILE__, __LINE__)
#define GetAlignResult()	__GetAlignResult(__FILE__, __LINE__)
#else
#define T(x)	\
	x *Get##x() \
		{ return (x *) GetObj(OT_##x); }
#endif
#include "objtypes.h"

	Obj *AllocNew(ObjType Type);
	void FreeObj(Obj *obj);
	uint GetFreeCount(ObjType Type) const;
	uint GetBusyCount(ObjType Type) const;
	uint GetMaxRefCount(ObjType Type) const;
	float GetTotalMem(ObjType Type) const;
	void LogBusy() const;

public:
	static void LogGlobalStats();
	static void LogGlobalBusy();

#if	DEBUG
	void ValidateType(ObjType Type) const;
	void Validate() const;
#endif
	};

#endif // objmgr_h
