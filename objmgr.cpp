#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "cpplock.h"

#undef Up
#undef Down

ObjMgr **ObjMgr::m_OMs;
unsigned ObjMgr::m_ThreadCount;

ObjMgr *ObjMgr::GetObjMgr()
	{
	LOCK();
	unsigned ThreadIndex = GetThreadIndex();
	if (ThreadIndex >= m_ThreadCount)
		{
		unsigned NewThreadCount = ThreadIndex + 32;
		ObjMgr **NewOMs = myalloc(ObjMgr *, NewThreadCount);
		zero_array(NewOMs, NewThreadCount);
		if (m_ThreadCount > 0)
			memcpy(NewOMs, m_OMs, m_ThreadCount*sizeof(ObjMgr *));
		m_OMs = NewOMs;
		m_ThreadCount = NewThreadCount;
		}
	if (m_OMs[ThreadIndex] == 0)
		m_OMs[ThreadIndex] = new ObjMgr(ThreadIndex);
	UNLOCK();
	return m_OMs[ThreadIndex];
	}

const char *ObjTypeToStr(ObjType Type)
	{
	switch (Type)
		{
#define T(x)	case OT_##x: return #x;
#include "objtypes.h"
		}
	return "OT_??";
	}

const char *ObjTypeToStr2(ObjType Type)
	{
	switch (Type)
		{
	case OT_SeqInfo:
		return "SI";
	case OT_PathInfo:
		return "PI";
	case OT_AlignResult:
		return "AR";
		}
	return "??";
	}

ObjMgr::ObjMgr(uint ThreadIndex)
	{
	m_ThreadIndex = ThreadIndex;
	zero_array(m_Free, OTCount);
	zero_array(m_Busy, OTCount);

#if	DEBUG || TRACE_OBJS
	m_Validate = false;
	zero_array(m_BusyCounts, OTCount);
	zero_array(m_GetCallCounts, OTCount);
	zero_array(m_AllocCallCounts, OTCount);
	zero_array(m_FreeCallCounts, OTCount);
#endif
	}

void ObjMgr::Down(Obj *pObj)
	{
	ObjMgr *OM = GetObjMgr();
	OM->ThreadDown(pObj);
	}

void ObjMgr::Up(Obj *pObj)
	{
	ObjMgr *OM = GetObjMgr();
	OM->ThreadUp(pObj);
	}

Obj *ObjMgr::StaticGetObj(ObjType Type)
	{
	ObjMgr *OM = GetObjMgr();
	return OM->ThreadGetObj(Type);
	}

Obj *ObjMgr::AllocNew(ObjType Type)
	{
#if	DEBUG
	++(m_AllocCallCounts[Type]);
#endif
	Obj *pObj = 0;
	switch (Type)
		{

#define	T(x)	case OT_##x: pObj = new x; break;
#include "objtypes.h"

	default:
		assert(false);
		}
	pObj->m_ThreadIndex = m_ThreadIndex;
	return pObj;
	}

Obj *ObjMgr::ThreadGetObj(ObjType Type)
	{
#if	DEBUG
	++(m_GetCallCounts[Type]);
	assert(Type < OTCount);
#endif
	Obj *NewObj = 0;
	if (m_Free[Type] == 0)
		NewObj = AllocNew(Type);
	else
		{
		NewObj = m_Free[Type];

		assert(NewObj->m_RefCount == 0);
		m_Free[Type] = m_Free[Type]->m_Fwd;
		if (m_Free[Type])
			m_Free[Type]->m_Bwd = 0;
		}

	if (m_Busy[Type] != 0)
		{
		assert(m_Busy[Type]->m_Bwd == 0);
		m_Busy[Type]->m_Bwd = NewObj;
		}
	NewObj->m_Fwd = m_Busy[Type];
	m_Busy[Type] = NewObj;

	assert(NewObj != 0);
	NewObj->m_RefCount = 1;

#if	DEBUG || TRACE_OBJS
	++(m_BusyCounts[Type]);
#endif
#if DEBUG
	if (m_Validate)
		Validate();
#endif

	return NewObj;
	}

void ObjMgr::FreeObj(Obj *obj)
	{
	assert(obj->m_RefCount == 0);

	ObjType Type = obj->m_Type;
	assert(Type < OTCount);

#if	DEBUG
	++(m_FreeCallCounts[Type]);
#endif

	if (obj == m_Busy[Type])
		m_Busy[Type] = obj->m_Fwd;

	Obj *Prev = obj->m_Bwd;
	Obj *Next = obj->m_Fwd;

	if (Prev != 0)
		Prev->m_Fwd = Next;
	if (Next != 0)
		Next->m_Bwd = Prev;

	if (m_Free[Type] != 0)
		{
		assert(m_Free[Type]->m_Bwd == 0);
		m_Free[Type]->m_Bwd = obj;
		}
	obj->m_Fwd = m_Free[Type];
	obj->m_Bwd = 0;
	m_Free[Type] = obj;

#if	DEBUG
	if (m_Validate)
		Validate();
#endif
	}

unsigned ObjMgr::GetFreeCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned Length = 0;
	const Obj *obj = m_Free[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		++Length;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Length;
	}

unsigned ObjMgr::GetBusyCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned Length = 0;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		++Length;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Length;
	}

unsigned ObjMgr::GetMaxRefCount(ObjType Type) const
	{
	assert(Type < OTCount);
	unsigned MaxRefCount = 0;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		if (obj->m_RefCount > MaxRefCount)
			MaxRefCount = obj->m_RefCount;
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return MaxRefCount;
	}

float ObjMgr::GetTotalMem(ObjType Type) const
	{
	assert(Type < OTCount);
	float Total = 0.0f;
	const Obj *obj = m_Busy[Type];
	for (;;)
		{
		if (obj == 0)
			break;
		Total += obj->GetMemBytes();
		assert(obj->m_Type == Type);
		obj = obj->m_Fwd;
		if (obj == 0)
			break;
		}
	return Total;
	}

#if	DEBUG
void ObjMgr::ValidateType(ObjType Type) const
	{
	unsigned NA = m_AllocCallCounts[Type];
	unsigned NF = m_FreeCallCounts[Type];

	unsigned nb = 0;
	for (const Obj *obj = m_Busy[Type]; obj; obj = obj->m_Fwd)
		{
		++nb;
		assert(nb <= NA);
		assert(obj->m_Type == Type);
		assert(obj->m_RefCount > 0);
		if (obj->m_Bwd)
			assert(obj->m_Bwd->m_Fwd == obj);
		if (obj->m_Fwd)
			assert(obj->m_Fwd->m_Bwd == obj);
		}

	unsigned nf = 0;
	for (const Obj *obj = m_Free[Type]; obj; obj = obj->m_Fwd)
		{
		++nf;
		assert(nf <= NF);
		assert(obj->m_RefCount == 0);
		assert(obj->m_Type == Type);
		if (obj->m_Bwd)
			assert(obj->m_Bwd->m_Fwd == obj);
		if (obj->m_Fwd)
			assert(obj->m_Fwd->m_Bwd == obj);
		}
	assert(nb + nf == NA);
	assert(nb == m_BusyCounts[Type]);
	}

void ObjMgr::Validate() const
	{
	Die("Validate!");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		ValidateType(Type);
		}
	}
#endif // DEBUG

#if	TRACE_OBJS

static bool g_HdrDone = false;

static void LogHdr()
	{
	if (g_HdrDone)
		return;

	Log("\n");
	Log("Op        UID  Ty  RefCnt  BsyCnt  Src\n");
	Log("---  --------  --  ------  ------  ---\n");
	g_HdrDone = true;
	}

Obj *ObjMgr::ThreadGetObj(ObjType Type, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	Obj *pObj = ThreadGetObj(Type);
	LOCK();
	static unsigned UID;
	pObj->m_UID = UID++;
	UNLOCK();
	pObj->m_SourceFileName = FileName;
	pObj->m_SourceLineNr = LineNr;

//	Log(" type = %u", pObj->m_Type);
	asserta((unsigned) pObj->m_Type < OTCount);

	Log("Get  %8u", pObj->m_UID);
	Log("  %2.2s", ObjTypeToStr2(pObj->m_Type));
	Log("  %6u", pObj->m_RefCount);
	Log("  %6u", m_BusyCounts[pObj->m_Type]);
	Log("  %s(%u)", BaseName(FileName), LineNr);
	Log("\n");

	return pObj;
	}

Obj *ObjMgr::StaticGetObj(ObjType Type, const char *FileName, unsigned LineNr)
	{
	ObjMgr *OM = GetObjMgr();
	Obj *pObj = OM->ThreadGetObj(Type, FileName, LineNr);
	return pObj;
	}

void ObjMgr::ThreadUp(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	ThreadUp(pObj);
	Log("Up   %8u  %2.2s  %6u  %6u  %s(%u)\n",
		pObj->m_UID,
		ObjTypeToStr2(pObj->m_Type),
		pObj->GetRefCount(),
		m_BusyCounts[pObj->m_Type],
		BaseName(FileName),
		LineNr);
	}

void ObjMgr::Up(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	ObjMgr *OM = GetObjMgr();
	OM->ThreadUp(pObj, FileName, LineNr);
	}

void ObjMgr::ThreadDown(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	LogHdr();
	ThreadDown(pObj);
	Log("Dwn  %8u  %2.2s  %6u  %6u  %s(%u)\n",
		pObj->m_UID,
		ObjTypeToStr2(pObj->m_Type),
		pObj->GetRefCount(),
		m_BusyCounts[pObj->m_Type],
		BaseName(FileName),
		LineNr);
	}

void ObjMgr::Down(Obj *pObj, const char *FileName, unsigned LineNr)
	{
	ObjMgr *OM = GetObjMgr();
	OM->ThreadDown(pObj, FileName, LineNr);
	}

void ObjMgr::LogBusy() const
	{
	Log("\n");
	Log("Busy:\n");
	Log("     UID  Ty  RefCnt  BsyCnt  Src\n");
	Log("--------  --  ------  ------  ---\n");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		for (const Obj *pObj = m_Busy[i]; pObj; pObj = pObj->m_Fwd)
			Log("%8u  %2.2s  %6u  %s(%u)\n",
				pObj->m_UID,
				ObjTypeToStr2(pObj->m_Type),
				pObj->GetRefCount(),
				pObj->m_SourceFileName,
				pObj->m_SourceLineNr);
		}
	}

#endif // TRACE_OBJS

static unsigned g_FreeCount[OTCount];
static unsigned g_BusyCount[OTCount];
static unsigned g_MaxRefCount[OTCount];
static float g_Mem[OTCount];

#if DEBUG
static unsigned g_GetCallCount[OTCount];
static unsigned g_AllocCallCount[OTCount];
static unsigned g_FreeCallCount[OTCount];
#endif

void ObjMgr::UpdateGlobalStats()
	{
	for (unsigned ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		ObjMgr *OM = m_OMs[ThreadIndex];
		if (OM != 0)
			OM->ThreadUpdateGlobalStats();
		}
	}

void ObjMgr::ThreadUpdateGlobalStats()
	{
	LOCK();
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		g_FreeCount[Type] += GetFreeCount(Type);
		g_BusyCount[Type] += GetBusyCount(Type);
		g_MaxRefCount[Type] = max(g_MaxRefCount[Type], GetMaxRefCount(Type));
		g_Mem[Type] += GetTotalMem(Type);
#if	DEBUG
		g_GetCallCount[Type] += m_GetCallCounts[Type];
		g_AllocCallCount[Type] += m_AllocCallCounts[Type];
		g_FreeCallCount[Type] += m_FreeCallCounts[Type];
#endif
		}
	UNLOCK();
	}

void ObjMgr::LogGlobalStats()
	{
	Log("\n");
	Log("            Type        Busy        Free         Mem   MaxRefCnt        Gets      Allocs       Frees\n");
	Log("----------------  ----------  ----------  ----------  ----------  ----------  ----------  ----------\n");
	for (unsigned i = 0; i < OTCount; ++i)
		{
		ObjType Type = (ObjType) i;
		const char *Name = ObjTypeToStr(Type);

		Log("%16.16s", Name);
		Log("  %10.10s", IntToStr(g_BusyCount[Type]));
		Log("  %10.10s", IntToStr(g_FreeCount[Type]));
		Log("  %10.10s", MemBytesToStr(g_Mem[Type]));
		Log("  %10u", g_MaxRefCount[Type]);
#if	DEBUG
		Log("  %10.10s", IntToStr(g_GetCallCount[Type]));
		Log("  %10.10s", IntToStr(g_AllocCallCount[Type]));
		Log("  %10.10s", IntToStr(g_FreeCallCount[Type]));
#endif
		Log("\n");
		}
	LogThreadStats();
#if	TRACE_OBJS
//	LogBusy();
#endif
	}

uint ObjMgr::GetBusyCount(uint Type) const
	{
	asserta(Type < OTCount);
	uint n = 0;
	for (const Obj *p = m_Busy[Type]; p; p = p->m_Fwd)
		++n;
	return n;
	}

uint ObjMgr::GetFreeCount(uint Type) const
	{
	asserta(Type < OTCount);
	uint n = 0;
	for (const Obj *p = m_Free[Type]; p; p = p->m_Fwd)
		++n;
	return n;
	}

void ObjMgr::LogThreadStats()
	{
	Log("ObjMgr::LogThreadStats() %u threads\n",
	  m_ThreadCount);
	for (uint ThreadIndex = 0; ThreadIndex < m_ThreadCount; ++ThreadIndex)
		{
		const ObjMgr *OM = m_OMs[ThreadIndex];
		if (OM != 0)
			OM->LogStats();
		}
	}

void ObjMgr::LogStats() const
	{
	Log("ObjMgr::LogStats() this=%p\n", this);
#define T(x)	Log("  %7u busy  %7u free  %s\n", \
	GetBusyCount(OT_##x), GetFreeCount(OT_##x), #x);
#include "objtypes.h"	
	}

void ObjMgr::ThreadDownByIndex(uint ThreadIndex, Obj *pObj)
	{
	asserta(ThreadIndex < m_ThreadCount);
	ObjMgr *OM = m_OMs[ThreadIndex];
	asserta(OM != 0);
	OM->ThreadDown(pObj);
	}
