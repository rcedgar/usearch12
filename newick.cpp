#include "myutils.h"
#include "tree.h"

#define TRACE 0

/***
Tokens in Newick files are:
	( ) : , ;
	string
	'string'
	"string"
	[ comment ]

We can't safely distinguish between identifiers and floating point
numbers at the lexical level (because identifiers may be numeric,
or start with digits), so both edge lengths and identifiers are
returned as strings.
***/

char Tree::GetCharFailOnEof(FILE *f)
	{
	int c = GetChar(f);
	if (c == EOF)
		Die("Parsing tree, unexpected end-of-file");
	return (char) c;
	}

int Tree::GetChar(FILE *f)
	{
	int c = fgetc(f);
	if (c == EOF)
		return c;
	char Ch = (char) c;

	++m_NewickColNr;
	if (Ch == '\n')
		{
		++m_NewickLineNr;
		m_NewickColNr = 1;
		}
	return Ch;
	}

void Tree::SkipWhite(FILE *f)
	{
	for (;;)
		{
		int c = GetChar(f);
		if (c == EOF)
			return;
		if (!isspace(c))
			{
			ungetc(c, f);
			return;
			}
		}
	}

const char *NTTToStr(NEWICK_TOKEN_TYPE NTT)
	{
	switch (NTT)
		{
#define c(x)	case NTT_##x: return #x;
	c(Unknown)
	c(Lparen)
	c(Rparen)
	c(Colon)
	c(Comma)
	c(Semicolon)
	c(String)
	c(SingleQuotedString)
	c(DoubleQuotedString)
	c(Comment)
#undef c
		}
	return "??";
	}

NEWICK_TOKEN_TYPE Tree::GetToken(FILE *f, string &Token)
	{
	NEWICK_TOKEN_TYPE NTT = GetTokenLo(f, Token);
#if TRACE
	Log("  %u(%u) token %s '%s'\n",
	  m_NewickLineNr, m_NewickColNr, NTTToStr(NTT), Token.c_str());
#endif
	return NTT;
	}

NEWICK_TOKEN_TYPE Tree::GetTokenLo(FILE *f, string &Token)
	{
	Token.clear();
	if (feof(f))
		return NTT_EOF;

// Skip leading white space
	SkipWhite(f);

	char c = GetCharFailOnEof(f);

// In case a single-character token
	Token = c;

	unsigned uBytesCopied = 0;
	NEWICK_TOKEN_TYPE TT;
	switch (c)
		{
	case '(':
		return NTT_Lparen;

	case ')':
		return NTT_Rparen;

	case ':':
		return NTT_Colon;

	case ';':
		return NTT_Semicolon;

	case ',':
		return NTT_Comma;

	case '\'':
		TT = NTT_SingleQuotedString;
		c = GetCharFailOnEof(f);
		break;

	case '"':
		TT = NTT_DoubleQuotedString;
		c = GetCharFailOnEof(f);
		break;

	case '[':
		TT = NTT_Comment;
		break;

	default:
		TT = NTT_String;
		break;
		}

// Discard char already added
	Token.clear();
	for (;;)
		{
		if (TT != NTT_Comment)
			Token += c;
		int ic = GetChar(f);
		if (ic == EOF)
			return TT;
		c = (char) ic;

		switch (TT)
			{
		case NTT_String:
			if (0 != strchr("():;,", c))
				{
				ungetc(c, f);
				return NTT_String;
				}
			if (isspace(c))
				return NTT_String;
			break;

		case NTT_SingleQuotedString:
			if ('\'' == c)
				return NTT_String;
			break;

		case NTT_DoubleQuotedString:
			if ('"' == c)
				return NTT_String;
			break;

		case NTT_Comment:
			if (']' == c)
				return GetToken(f, Token);
			break;

		default:
			Die("Tree::GetToken, invalid TT=%u", TT);
			}
		}
	}

void Tree::FromNewickFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);
	m_NewickLineNr = 1;
	m_NewickColNr = 0;
	m_NewickDepth = 0;

	Clear();

	unsigned Root = GetGroupFromFile(f, UINT_MAX);
	asserta(Root == 0 || Root == THIRD_NODE);
	if (Root == THIRD_NODE)
		{
		uint ThirdRootChild = GetGroupFromFile(f, THIRD_NODE);
		string Token;
		NEWICK_TOKEN_TYPE NTT = GetToken(f, Token);
		if (NTT != NTT_Rparen)
			Die("(%u(%u) Expected ')' after third group, got '%s'",
			 m_NewickLineNr, m_NewickColNr, Token.c_str());
		
	//	uint FirstRootChild = m_Lefts[0];
		uint SecondRootChild = m_Rights[0];

	// Insert new node parent=root left=second_child right=third_child
	// This new node becomes second child of root.
	// first child of root unchanged.
		unsigned NewNodeIndex = SIZE(m_Parents);
		m_Parents.push_back(0);
		m_Labels.push_back("_unroot_");
		m_Lefts.push_back(SecondRootChild);
		m_Rights.push_back(ThirdRootChild);
		m_Lengths.push_back(DBL_MAX);

		m_Rights[0] = NewNodeIndex;
		m_Parents[SecondRootChild] = NewNodeIndex;
		m_Parents[ThirdRootChild] = NewNodeIndex;
		}

	string Token;
	NEWICK_TOKEN_TYPE NTT = GetToken(f, Token);

	if (NTT != NTT_Semicolon)
		Die("%u(%u) Expected ';' after root group, got '%s'",
		  m_NewickLineNr, m_NewickColNr, Token.c_str());

	CloseStdioFile(f);
	Validate();
	}

unsigned Tree::GetGroupFromFile(FILE *f, unsigned Parent)
	{
#if	TRACE
	Log("GetGroupFromFile parent=%u depth=%u\n", Parent, m_NewickDepth);
#endif

	++m_NewickDepth;
	string Token;
	NEWICK_TOKEN_TYPE NTT = GetToken(f, Token);

	unsigned NodeIndex = SIZE(m_Parents);
	m_Parents.push_back(Parent);
	m_Labels.push_back("");
	m_Lefts.push_back(UINT_MAX);
	m_Rights.push_back(UINT_MAX);
	m_Lengths.push_back(DBL_MAX);

// Group is either leaf name or (left, right).
	if (NTT == NTT_String)
		m_Labels[NodeIndex] = Token;
	else if (NTT == NTT_Lparen)
		{
#if TRACE
		Log("Get Left(%u, '%s')\n", NodeIndex, m_Labels[NodeIndex].c_str());
#endif
		unsigned Left = GetGroupFromFile(f, NodeIndex);
#if TRACE
		Log("  Left(%u, '%s')=%u, '%s'\n", NodeIndex, m_Labels[NodeIndex].c_str(), Left, m_Labels[Left].c_str());
#endif
		NTT = GetToken(f, Token);
		if (NTT != NTT_Comma)
			Die("Newick file line %u col %u, expected ',', got '%s'",
			  m_NewickLineNr, m_NewickColNr, Token.c_str());
#if TRACE
		Log("Get Right(%u, '%s')\n", NodeIndex, m_Labels[NodeIndex].c_str());
#endif
		unsigned Right = GetGroupFromFile(f, NodeIndex);
#if TRACE
		Log("  Right(%u, '%s')=%u, '%s'\n", NodeIndex, m_Labels[NodeIndex].c_str(), Right, m_Labels[Right].c_str());
#endif

		m_Lefts[NodeIndex] = Left;
		m_Rights[NodeIndex] = Right;

	// Expect closing parenthesis (or comma if > 2-ary)
		NTT = GetToken(f, Token);
		if (NTT == NTT_Comma)
			{
			// Die("Tree is not binary");
			return THIRD_NODE;
			}
		else if (NTT != NTT_Rparen)
			Die("Newick file line %u col %u, expected ')', got '%s'",
			  m_NewickLineNr, m_NewickColNr, Token.c_str());
		}
	else
		Die("Newick file line %u col %u, expected '(' or leaf name, got '%s'",
		  m_NewickLineNr, m_NewickColNr, Token.c_str());
	++m_NewickDepth;

// Group optionally followed by colon & edge length or
// internal node name then colon and edge length
	SkipWhite(f);
	if (feof(f))
		return NodeIndex;

	char c = GetCharFailOnEof(f);
	if (isalnum(c) || c == '\'')
		{
		asserta(m_Labels[NodeIndex] == "");
		ungetc(c, f);
		NTT = GetToken(f, Token);
		if (NTT != NTT_String)
			Die("%u(%u), expected internal node name, got '%s'", 
			  m_NewickLineNr, m_NewickColNr, Token.c_str());
		m_Labels[NodeIndex] = Token;
		c = GetCharFailOnEof(f);
		}

	if (c == ':')
		{
		NTT = GetToken(f, Token);
		if (NTT != NTT_String)
			Die("%u(%u), expected edge length, got '%s'", 
			  m_NewickLineNr, m_NewickColNr, Token.c_str());
		double EdgeLength = StrToFloat(Token);
		m_Lengths[NodeIndex] = EdgeLength;
		return NodeIndex;
		}

	ungetc(c, f);
	return NodeIndex;
	}
