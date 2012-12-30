#ifndef BASIC
#define BASIC

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <map>


#define IN
#define OUT
#define INVALID -1

using namespace std;

template<class T> class FSTVector:public std::vector<T>
{
public:
	void push_back(const T &obj)
	{
		if(std::vector<T>::capacity() <= 10 + std::vector<T>::size())
		{
			std::vector<T>::reserve(int(1+std::vector<T>::capacity() * 1.1)); // 减小内存使用
		}
		std::vector<T>::push_back(obj);
	}
};

struct FSTARC
{
	int state_i;
	int state_o;
	float weight; // 既可以是P，也可以是B

	string input;
	string output;
	FSTARC() {}

	FSTARC(int statein, int stateout, string i, string o = "", float weight = 0)
	{
		state_i = statein;
		state_o = stateout;
		this->weight = weight;
		input = i;
		output = o;
	}

	void set(int statein, int stateout, string io, float weight = 0)
	{
		state_i = statein;
		state_o = stateout;
		this->weight = weight;
		input = io;
		output = io;
	}

	bool isValid(){
		return (state_i!=INVALID || state_o!=INVALID);
	}
};

struct FSTSTATE
{
	//  FSTVector<int> arc_i;
	FSTVector<int> arc_o; // 约定第0个位置存放回退的弧
	FSTSTATE() {}
	FSTSTATE(OUT int o)
	{
		arc_o.push_back(INVALID);//用-1占用0号输出弧，代表无回退弧
		arc_o.push_back(o);
	}
	bool isValid(){
		return !arc_o.empty();
	}
};

typedef FSTVector<string> vs;


class FST
{
	/*建立FST相关的函数以及FST中的节点变量*/
public:
	int state_s; // 开始状态
	int state_e; // 结束状态, 如果FST由Ngram构建，则state_s=state_e=0（空状态是整个Ngram的开始，也是结束）
	vector<vector <int> > arcs_in; //用于暂存入弧，用完即删
	FSTVector<FSTSTATE> states; // 约定第0个状态是空状态
	FSTVector<FSTARC> arcs;

public:
	FST();
	~FST();

public:
	friend bool operator == (const FSTARC &arc1, const FSTARC &arc2);
	friend vector<int> &operator += (vector<int> &dst , const vector<int> &src);
	bool ConvertNgram2FST(IN string filename);  // 请调用 AddUniGram,AddBiGram,AddTriGram 实现
	bool AddUniGram(float prob, string word, float backoff );
	bool AddBiGram(float prob, string word1, string word2, float backoff);
	bool AddTriGram(float prob, string word1, string word2, string word3, float backoff);
	float GetSentenceScore(IN const FST &fst, IN string sentence);
	bool WriteFST(string filename);
	bool FindWord(int from_state, string word, int &next_state, int &next_arc);
	void UnionFromFST(const FST &fst1, const FST &fst2);
	FST &operator + (const FST &fst) {
		UnionFromFST(*this, fst);
		return *this;
	}
	void ShowInfo() {
		printf("sizeof(FSTSTATE)=%d count=%d/%d\n", (int)sizeof(FSTSTATE), (int)states.size(), (int)states.capacity());
		printf("sizeof(FSTARC)=%d count=%d/%d\n", (int)sizeof(FSTARC), (int)arcs.size(), (int)arcs.capacity());
	}
	static void split(string line, vs &parts, string pattern);
	void slim();
	float ConvertInputToOutput(string input, int beamwidth, string & output);
	float Transduce(int cur_state, string phone, int max_arc_num, int * arc_ids, float * arc_weights, int & return_arc_num);
	void  UpdateArcsI();
	
};

struct PATH
{
	float score;
	int state_id;
	int arc_id;
	string output;

	bool isValid()
	{
		return !(state_id==INVALID);
	}

	PATH()
	{
		score=0;
		state_id=INVALID;//state_id=-1时代表无效状态
		arc_id=0;
		output="";
	}

	void AddOutput(const string & word)//输出到beam上
	{
		this->output+=word;
	}
};

map<string, vector<string>	> confusion;

/////////////////////////////functions/////////////////////////////////////////////////////////

FST::FST(){}

FST::~FST(){
	states.clear();
	arcs.clear();
}

bool FST::WriteFST(string filename)
{
	fstream fst(filename.c_str(), ios::out);
	if (!fst)
	{
		return false;
	}
	// fst << "#FSTBasic MaxPlus" << endl;
	// fst << "I 0" << endl;
	// fst << "F 0" << endl;

	for (int i = 0; i < arcs.size(); i++)
	{
		FSTARC &arc = arcs[i];
		string input = arc.input.size() > 0 ? arc.input : ",";
		string output = arc.output.size() > 0 ? arc.output : ",";
		// fst << "T " << arc.state_i << " " << arc.state_o << " " << input << " " << output << " " << arc.weight << endl;
		fst << arc.state_i << " " << arc.state_o << " " << input << " " << output << " " << arc.weight << endl;
		//      fst<<"T "<<arc.state_o<<" "<<input<<" "<<output<<" "<<arc.weight<<endl;
	}
	fst << "0 1" << endl; //added later
	fst.close();
	return true;
}


/*读入Ngrams，转换成WFST*/
bool FST::ConvertNgram2FST(IN string filename) // 请调用 AddUniGram,AddBiGram,AddTriGram 实现
{

	fstream dictfile;
	string line;
	vs substrs;
	float P, B;
	int linenum = 0;
	string tmpstr, word, word1, word2, word3;
	state_s = 0;
	state_e = 0;
	dictfile.open(filename.c_str(), ios::in);
	if (!dictfile) {
		cout << "Sorry,can't open file." << endl;
		return 1;
	}

	FSTSTATE headstate;
	headstate.arc_o.push_back(INVALID);
	states.push_back(headstate);//默认创建一个ε头节点,无回退弧（以-1表示）

	while (!dictfile.eof())
	{
		getline(dictfile, line);
		if (!line.compare("\\data\\")) break;
		//寻找"\data\"出现
	}
	getline(dictfile, line);
	split(line, substrs, "=");
	int n_1gram = atoi(substrs[1].c_str());
	getline(dictfile, line);
	split(line, substrs, "=");
	int n_2gram = atoi(substrs[1].c_str());
	getline(dictfile, line);
	split(line, substrs, "=");
	int n_3gram = atoi(substrs[1].c_str());

	/////1-gram////
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		//      line.resize(8);
		if (!line.compare("\\1-grams:")) break;
		//寻找"\1-grams"出现
	}
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		linenum++; if (linenum % 1000 == 0) printf("%d//%d: %s\n", linenum, n_1gram, line.c_str());
		if (!line.compare("")) break;
		//处理1-grams
		split(line, substrs, "\r\n\t ");
		tmpstr = substrs[0];
		P = (float)atof(tmpstr.c_str());
		word = substrs[1];
		if (substrs.size() > 2) //初始化B
		{
			tmpstr = substrs[2];
			B = (float)atof(tmpstr.c_str());
		}
		else
		{
			B = 0; //B默认为1
		}
		AddUniGram(P, word, B);
	}
	////2-gram/////
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		//      line.resize(8);
		if (!line.compare("\\2-grams:")) break;
		//寻找"\2-grams"
	}
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		linenum++; if (linenum % 1000 == 0) printf("%d//%d: %s\n", linenum, n_2gram, line.c_str());
		if (!line.compare("")) break;
		//处理2-grams
		split(line, substrs, "\r\n\t ");
		tmpstr = substrs[0];
		P = (float)atof(tmpstr.c_str());
		word1 = substrs[1];
		word2 = substrs[2];
		if (substrs.size() > 3)
		{
			tmpstr = substrs[3];
			B = (float)atof(tmpstr.c_str());
		}
		else
		{
			B = 0;
		}
		AddBiGram(P, word1, word2, B);
	}
	////3-gram////
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		//      line.resize(8);
		if (!line.compare("\\3-grams:")) break;
		//寻找"\3-grams"
	}
	while (!dictfile.eof())
	{
		getline(dictfile, line);
		linenum++; if (linenum % 1000 == 0) printf("%d//%d: %s\n", linenum, n_3gram, line.c_str());
		if (!line.compare("")) break;
		//处理3-grams
		split(line, substrs, "\r\n\t ");
		tmpstr = substrs[0];
		P = (float)atof(tmpstr.c_str());
		word1 = substrs[1];
		word2 = substrs[2];
		word3 = substrs[3];
		if (substrs.size() > 4)
		{
			tmpstr = substrs[4];
			B = (float)atof(tmpstr.c_str());
		}
		else
		{
			B = 0;
		}
		AddTriGram(P, word1, word2, word3, B);
	}
	return 0;
}


// 1.1 添加1元语言模型到WFST中，注意如果WFST为空应该如何添加UniGram
bool FST::AddUniGram(float prob, string word, float backoff )
{
	FSTARC tmparc;
	FSTSTATE tmpstate;

	tmparc.set(states.size(), 0, "", backoff);

	arcs.push_back(tmparc);
	tmpstate.arc_o.push_back(arcs.size() - 1); //创建回退弧

	tmparc.set(0, states.size(), word, prob);
	//  tmpstate.arc_i.push_back(arcs.size());
	states[0].arc_o.push_back(arcs.size());
	arcs.push_back(tmparc);//创建前进弧

	//tmpstate.state_word=word;
	states.push_back(tmpstate);//插入新节点
	return 1;
}

// 1.2 添加2元语言模型到WFST中，注意如果word1的UniGram没有找到，返回false
bool FST::AddBiGram(float prob, string word1, string word2, float backoff)
{
	FSTVector<FSTSTATE>::iterator it_state;
	static int UniGram_counter = states.size() - 1; //UniGram的数量为当前的总states数-1;
	int state_back, state_prev, arc_id;

	if (!FindWord(0, word1, state_prev, arc_id))
		return false;
	if (!FindWord(0, word2, state_back, arc_id))
		return false;

	FSTARC tmparc;
	FSTSTATE tmpstate;
	////////////////往后添加BiGram FST//////////
	tmparc.set(states.size(), state_back, "", backoff);
	arcs.push_back(tmparc);
	tmpstate.arc_o.push_back(arcs.size() - 1); //创建回退弧

	tmparc.set(state_prev, states.size(), word2, prob);
	//tmpstate.arc_i.push_back(arcs.size());
	states[state_prev].arc_o.push_back(arcs.size());
	arcs.push_back(tmparc);//创建前进弧

	//tmpstate.state_word=word1+word2;
	states.push_back(tmpstate);//插入新节点
	return true;
}

// 1.3 添加3元语言模型到WFST中，注意如果(word1,word2)的BiGram没有找到，应该先补上
bool FST::AddTriGram(float prob, string word1, string word2, string word3, float backoff)
{
	FSTVector<FSTSTATE>::iterator it_state;
	int stateid, arcid, arc_prev, arc_back;
	int state_back, state_prev; //存储回退节点ID和上一个节点ID；
	string word_prev = word1 + word2;
	string word_back = word3; //存储回退节点的状态

	if (!FindWord(0, word1, stateid, arcid))
		return false;

	if (!FindWord(stateid, word2, state_prev, arc_prev))
	{
		int arc_back = states[stateid].arc_o[0];

		float P = arcs[arc_back].weight;

		if (!FindWord(0, word2, state_prev, arc_prev))
			return false;

		P += arcs[arcid].weight;

		AddBiGram(P, word1, word2, 0.0f);

		state_prev = states.size() - 1; //记录新建立的BiGram的节点ID
	}

	if (!FindWord(0, word2, stateid, arcid))
		return false;

	if (!FindWord(stateid, word3, state_back, arc_back))
	{
		if (!FindWord(0, word3, state_back, arc_back))
			return false;
	}

	//此时回退节点和上一个节点的ID已经确定

	FSTARC tmparc;
	FSTSTATE tmpstate;

	tmparc.set(states.size(), state_back, "", backoff);
	arcs.push_back(tmparc);
	tmpstate.arc_o.push_back(arcs.size() - 1); //创建回退弧

	tmparc.set(state_prev, states.size(), word3, prob);
	//tmpstate.arc_i.push_back(arcs.size());
	states[state_prev].arc_o.push_back(arcs.size());
	arcs.push_back(tmparc);//创建前进弧

	// tmpstate.state_word=word1+word2+word3;
	states.push_back(tmpstate);//插入新节点
	return true;

	return 1;
}

void FST::split(string line, vs &parts, string pattern)
{
	int n = line.size() + 1;
	char *buffer = new char[n];
	strcpy(buffer, line.c_str());
	char *p = (char *)buffer, *q = 0 ;

	const char *seg = "\r \n\t";
	if (pattern.size() > 0)
		seg = pattern.c_str();

	parts.clear();

	while (1 == 1)
	{
		// trim left
		while (*p && strchr(seg, *p))p++;
		if (*p == 0) break;

		// get info
		q = p;
		while (*q && !strchr(seg, *q))q++;
		char ch = *q; * q = 0;
		parts.push_back(p);
		*q = ch; p = q;
	}

	delete [] buffer;
}

bool FST::FindWord(int from_state, string word, int &next_state, int &next_arc)
{
	//  if (next_state==INVALID || next_arc==INVALID)
	//  {
	//      return false;
	//  }


	int st = 1; // state 0 is eps state
	int ed = states[from_state].arc_o.size() - 1;
	int md = (st + ed) / 2;

	next_state = INVALID;
	next_arc = INVALID;
	while (st <= ed)
	{
		next_arc = states[from_state].arc_o[md];
		int res_md = arcs[next_arc].input.compare(word);
		if (res_md == 0)
		{
			next_state = arcs[next_arc].state_o;
			break;
		}

		if (res_md < 0)
		{
			st = md + 1;
			md = (st + ed) / 2;
		}
		else
		{
			ed = md - 1;
			md = (st + ed) / 2;
		}
	}

	if (next_state < 0)
	{
		// printf("WARNING: cannot find forward %s\n", word.c_str());

		st = 0;
		ed = states[from_state].arc_o.size() - 1;
		for (md = ed; md >= st;  md --)
		{
			next_arc = states[from_state].arc_o[md];
			if (!arcs[next_arc].input.compare(word))
			{
				next_state = arcs[next_arc].state_o;//获得所对应的STATE位置(int)state_back
				break;
			}
		}
	}

	if (next_state < 0)
	{
		return false;
	}

	return true;
}

bool operator == (const FSTARC &arc1, const FSTARC &arc2)
{
	if (arc1.input == arc2.input && arc1.output == arc2.output && arc1.weight == arc2.weight)
	{
		return true;
	}
	else
	{
		return false;
	}
}

vector<int> operator + (const vector<int> &dst, const  vector<int> &src)
{
	vector<int> tmp = dst;
	for (int i = 0; i < src.size(); i++)
	{
		tmp.push_back(src[i]);
	}
	return tmp;
}

vector<int> &operator += (vector<int> &dst , const vector<int> &src)
{
	for (int i = 0; i < src.size(); i++)
	{
		dst.push_back(src[i]);
	}
	return dst;
}


void FST::slim()
{

	// 1. 删除多余的弧
	int src_arc_id, dst_arc_id;
	for (src_arc_id = dst_arc_id = 0; src_arc_id < arcs.size(); src_arc_id++)
	{
		if (arcs[src_arc_id].isValid())
		{
			arcs[dst_arc_id ++] = arcs[src_arc_id];
		}
	}

	arcs.resize(dst_arc_id);

	// 2. 重建所有状态的输出弧
	bool NullBackoffState;
	for (int i = 0; i < states.size(); i++)
	{
		if (states[i].arc_o.empty())
		{
			continue;//无效状态依然保留
		}
		NullBackoffState = (states[i].arc_o[0] == INVALID);
		states[i].arc_o.clear();
		if (NullBackoffState)
		{
			states[i].arc_o.push_back(INVALID);
		}
	}
	//  for(int i=0;i<states.size();i++)
	//  {
	//      states[i].arc_o.clear();
	//  }
	for (int i = 0; i < arcs.size(); i++)
	{
		int state_id = arcs[i].state_i;
		states[state_id].arc_o.push_back(i);
	}

	//3.删除多余的状态
	UpdateArcsI();
	int src_state_id, dst_state_id;
	for (src_state_id = dst_state_id = 0; src_state_id < states.size(); src_state_id++)
	{
		if (states[src_state_id].isValid())
		{
			states[dst_state_id] = states[src_state_id];
			arcs_in[dst_state_id++] = arcs_in[src_state_id];
		}
	}
	states.resize(dst_state_id);
	arcs_in.resize(dst_state_id);


	//4.修改弧的出入状态
	for (int i = 0; i < states.size(); i++)
	{
		vector<int> &arc_outs = states[i].arc_o;
		vector<int> &arc_ins = arcs_in[i];
		if (arc_outs[0] != INVALID)
		{
			arcs[arc_outs[0]].state_i = i;
		}
		for (int j = 1; j < arc_outs.size(); j++)
		{
			arcs[arc_outs[j]].state_i = i;
		}

		for (int j = 0; j < arc_ins.size(); j++)
		{
			arcs[arc_ins[j]].state_o = i;
		}
	}
}

void FST::UpdateArcsI()
{
	arcs_in.resize(states.size());
	for (int i = 0; i < arcs.size(); i++)
	{
		if (arcs[i].state_o != INVALID)
			arcs_in[arcs[i].state_o].push_back(i);
	}
}


#endif
