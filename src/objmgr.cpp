#include "myutils.h"
#include "objmgr.h"
#include "seqinfo.h"
#include "pathinfo.h"
#include "alignresult.h"
#include "mymutex.h"
#include "obj.h"

void Obj::Up()
	{
	m_Owner->Up(this);
	}

void Obj::Down()
	{
	m_Owner->Down(this);
	}

vector<ObjMgr *> ObjMgr::m_OMs;
vector<ObjMgr *> ObjMgr::m_FreeOMs;

void ObjMgr::FreeObjMgr(ObjMgr *OM)
	{
	static MUTEX(mut, "ObjMgr::FreeObjMgr");
	mut.lock();
	m_FreeOMs.push_back(OM);
	vector<ObjMgr *> NewOMs;
	for (size_t i = 0; i < m_OMs.size(); ++i)
		{
		ObjMgr *OM2 = m_OMs[i];
		if (OM2 != OM)
			NewOMs.push_back(OM2);
		}
	m_OMs = NewOMs;
	mut.unlock();
	}

ObjMgr *ObjMgr::CreateObjMgr()
	{
	static MUTEX(mut, "ObjMgr::CreateObjMgr");
	mut.lock();
	ObjMgr *OM = 0;
	if (!m_FreeOMs.empty())
		{
		OM = m_FreeOMs.back();
		m_FreeOMs.pop_back();
		}
	else
		OM = new ObjMgr;
	m_OMs.push_back(OM);
	mut.unlock();
	return OM;
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

ObjMgr::ObjMgr()
	{
#if DEBUG
	m_Validate = false;
#endif
	zero_array(m_Free, OTCount);
	zero_array(m_Busy, OTCount);
	}

void ObjMgr::Up(Obj *pObj)
	{
#if	DEBUG
	if (m_Validate)
		Validate();
#endif

	assert(pObj->m_RefCount != 0);
	++pObj->m_RefCount;

#if	DEBUG
	if (m_Validate)
		Validate();
#endif
	}

void ObjMgr::Down(Obj *pObj)
	{
#if	DEBUG
	if (m_Validate)
		Validate();
#endif
	assert(pObj->m_RefCount > 0);
	--pObj->m_RefCount;
	if (pObj->m_RefCount == 0)
		{
		ObjType Type = pObj->m_Type;
		FreeObj(pObj);
		pObj->OnZeroRefCount();
		}

#if	DEBUG
	if (m_Validate)
		Validate();
#endif
	}

Obj *ObjMgr::AllocNew(ObjType Type)
	{
	Obj *pObj = 0;
	switch (Type)
		{

#define	T(x)	case OT_##x: pObj = new x; break;
#include "objtypes.h"

	default:
		assert(false);
		}
	pObj->m_Owner = this;
	return pObj;
	}

#if TRACK_OBJS
Obj *ObjMgr::GetObj(ObjType Type, const char *SrcFN, uint LineNr)
#else
Obj *ObjMgr::GetObj(ObjType Type)
#endif
	{
	assert(Type < OTCount);
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

#if DEBUG
	if (m_Validate)
		Validate();
#endif
#if TRACK_OBJS
	NewObj->m_SourceFileName = SrcFN;
	NewObj->m_SourceLineNr = LineNr;
#endif

	return NewObj;
	}

void ObjMgr::FreeObj(Obj *obj)
	{
	assert(obj->m_RefCount == 0);

	ObjType Type = obj->m_Type;
	assert(Type < OTCount);

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
	unsigned nb = 0;
	for (const Obj *obj = m_Busy[Type]; obj; obj = obj->m_Fwd)
		{
		++nb;
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
		assert(obj->m_RefCount == 0);
		assert(obj->m_Type == Type);
		if (obj->m_Bwd)
			assert(obj->m_Bwd->m_Fwd == obj);
		if (obj->m_Fwd)
			assert(obj->m_Fwd->m_Bwd == obj);
		}
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

void Obj *pObj, const char *FileName, unsigned LineNr->Up()
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

void Obj *pObj, const char *FileName, unsigned LineNr->Down()
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

void ObjMgr::LogBusy() const
	{
#if TRACK_OBJS
	for (uint iType = 0; iType < OTCount; ++iType)
		{
		for (const Obj *pObj = m_Busy[iType]; pObj; pObj = pObj->m_Fwd)
			Log("%s:%d\n", pObj->m_SourceFileName, pObj->m_SourceLineNr);
		}
#else
	Log("ObjMgr::LogBusy(), !TRACK_OBJS\n");
#endif
	}

void ObjMgr::LogStats() const
	{
	Log("ObjMgr::LogStats() this=%p\n", this);
#define T(x)	Log("  %7u busy  %7u free  %s\n", \
	GetBusyCount(OT_##x), GetFreeCount(OT_##x), #x);
#include "objtypes.h"	
	}

void ObjMgr::LogGlobalBusy()
	{
#if TRACK_OBJS
	uint N = SIZE(m_OMs);
	for (uint i = 0; i < N; ++i)
		m_OMs[i]->LogBusy();
#endif
	}

void ObjMgr::LogGlobalStats()
	{
	uint N = SIZE(m_OMs);
	Log("ObjMgr::LogGlobalStats(), %u OMs\n", N);
	Log("");
	Log("     Obj        Busy        Free   MaxRefCnt         Mem\n");
	Log("--------  ----------  ----------  ----------  ----------\n");
	for (uint iType = 0; iType < OTCount; ++iType)
		{
		ObjType Type = (ObjType) iType;
		uint FreeCount = 0;
		uint BusyCount = 0;
		uint MaxRefCount = 0;
		float TotalMem = 0;
		for (uint i = 0; i < N; ++i)
			{
			FreeCount += m_OMs[i]->GetFreeCount(Type);
			BusyCount += m_OMs[i]->GetBusyCount(Type);
			TotalMem += m_OMs[i]->GetTotalMem(Type);
			MaxRefCount = max(MaxRefCount, m_OMs[i]->GetMaxRefCount(Type));
			}
		Log("%8.8s  %10u  %10u  %10u  %10.10s\n",
			ObjTypeToStr(Type), BusyCount, FreeCount,
			MaxRefCount, MemBytesToStr(TotalMem));
		}
#if TRACK_OBJS
	LogGlobalBusy();
#endif
	}
