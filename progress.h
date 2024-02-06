#pragma once

class SeqSource;

void StopProgressThread();
void StartProgressThread();

typedef void (PTR_PROGRESS_CB)(string &s);

void ProgressNote(const char *fmt, ...);
void ProgressNoteLog(const char *fmt, ...);
void ProgressStart(const char *fmt, ...);
void ProgressStartCB(PTR_PROGRESS_CB CB, const char *fmt, ...);
void ProgressStartSS(SeqSource &SS, const char *fmt, ...);
void ProgressDone();
void ProgressLoop(uint64 *ptri, uint64 N, const char *fmt, ...);
void ProgressLoop(uint32 *ptri, uint32 N, const char *fmt, ...);
