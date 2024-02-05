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

private:
	ObjMgr();

#if	DEBUG || TRACE_OBJS
	bool m_Validate;
	uint m_BusyCounts[OTCount];
	uint m_GetCallCounts[OTCount];
	uint m_FreeCallCounts[OTCount];
	uint m_AllocCallCounts[OTCount];
#endif

public:
	static vector<ObjMgr *> m_OMs;

public:
	static ObjMgr *CreateObjMgr();

public:
	void Up(Obj *pObj);
	void Down(Obj *pObj);
	uint GetBusyCount(uint Type) const;
	uint GetFreeCount(uint Type) const;
	void LogStats() const;
	Obj *GetObj(ObjType Type);
	//Obj *GetObj(ObjType Type, const char *FileName, uint LineNr);

#define T(x)	\
	x *Get##x() \
		{ return (x *) GetObj(OT_##x); } \
	//x *Get##x(const char *FileName, uint LineNr) \
	//	{ return (x *) GetObj(OT_##x, FileName, LineNr); }
#include "objtypes.h"

	Obj *AllocNew(ObjType Type);
	void FreeObj(Obj *obj);
	uint GetFreeCount(ObjType Type) const;
	uint GetBusyCount(ObjType Type) const;
	uint GetMaxRefCount(ObjType Type) const;
	float GetTotalMem(ObjType Type) const;

public:
	static void LogGlobalStats();

#if	DEBUG
	void ValidateType(ObjType Type) const;
	void Validate() const;
#endif
	};

#endif // objmgr_h
