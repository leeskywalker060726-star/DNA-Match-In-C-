#include<bits/stdc++.h>
#include<chrono>
using namespace std;
using namespace std::chrono;
#define N 30005

//FRec: Fasta序列记录(Fasta Record)
struct FRec{
    string h;//h:header
    string s;//s:sequence
};

//Res: 存储比对结果(Result)
struct Res{
    int qs,qe,rs,re;
};

//rdFasta: 读取fasta文件(read Fasta)
vector<FRec>rdFasta(const string&fn){ //fn:filename
    vector<FRec>v;
    ifstream f(fn);
    string l,h="",s=""; //l:line, h:header, s:seq
    while(getline(f,l)){
        if(!l.empty()&&l.back()=='\r'){
            l.pop_back();
        }
        if(l.empty()){
            continue;
        }
        if(l[0]=='>'){
            if(!h.empty()){
                v.push_back({h,s});
                s="";
            }
            h=l.substr(1);
        }else{
            s+=l;
        }
    }
    if(!h.empty()){
        v.push_back({h,s});
    }
    f.close();
    return v;
}

//getRC: 获取反向互补序列(get Reverse Complement)
string getRC(const string&s){ //s:seq
    string r=""; //r:rc seq
    for(int i=(int)s.length()-1;i>=0;--i){
        switch(s[i]){
            case 'A':{r+='T';break;}
            case 'T':{r+='A';break;}
            case 'C':{r+='G';break;}
            case 'G':{r+='C';break;}
            default:{r+=s[i];break;}
        }
    }
    return r;
}

string R,Q; //R:Ref, Q:Query
int f[N],g[N],lR,lQ; //lR:lenR, lQ:lenQ. f:dp分值, g:dp回溯
int grs[N]; //grs:g_r_start, 记录回溯对应的R起始位置

//getPM: 获取完美匹配位置(get Perfect Match)
int getPM(int l,int r){
    if(r-l+1>lR){
        return-1;
    }
    string h=Q.substr(l-1,r-l+1);
    int p1=R.find(h); //p1:pos1
    if(p1!=string::npos){
        return p1;
    }
    string h2=getRC(h);
    int p2=R.find(h2); //p2:pos2
    if(p2!=string::npos){
        return p2;
    }
    return-1;
}

int main(){
    string fn="Lab1.fasta";
    vector<FRec>d=rdFasta(fn); //d:data
    if(d.size()<4){
        return 0;
    }
    
    Q=d[2].s; 
    R=d[3].s;
    lR=(int)R.size();
    lQ=(int)Q.size();
    
    for(int i=1;i<=lQ;++i){
        f[i]=i;
    }
    f[0]=0;
    
    for(int i=0;i<=lQ;++i){
        f[i]=0; 
        grs[i]=-1; 
    }
    
    //fw: forward匹配数, rv: reverse匹配数
    static int fw[N];
    static int rv[N];
    auto st=high_resolution_clock::now(); //st:start time

    for(int i=0;i<lQ;++i){
        if(f[i]>f[i+1]){
            f[i+1]=f[i];
            g[i+1]=i; 
            grs[i+1]=-1; 
        }

        memset(fw,0,sizeof(int)*(lR+1));
        memset(rv,0,sizeof(int)*(lR+1));

        for(int j=1;j<=min(lQ-i,240);++j){
            char fc=Q[i+j-1]; //fc:fwd char
            char rc; //rc:rev char
            switch(fc){
                case 'A':{rc='T';break;}
                case 'T':{rc='A';break;}
                case 'C':{rc='G';break;}
                case 'G':{rc='C';break;}
                default:{rc=fc;break;}
            }

            int mF=1e9,pF=-1; //mF:min fwd edits, pF:pos fwd
            int mR=1e9,pR=-1; //mR:min rev edits, pR:pos rev
            
            for(int k=0;k<=lR-j;++k){
                if(R[k+j-1]==fc){
                    fw[k]++;
                }
                rv[k]=rv[k+1]+(R[k]==rc);

                if(j>=30){
                    int eF=j-fw[k]; //eF:edits fwd
                    if(eF<mF){
                        mF=eF;
                        pF=k;
                    }
                    
                    int eR=j-rv[k]; //eR:edits rev
                    if(eR<mR){
                        mR=eR;
                        pR=k;
                    }
                }
            }

            if(j>=30){
                int ed=1e9; //ed:final edits
                int rp=-1; //rp:r start pos

                if(mF==0){
                    ed=0;
                    rp=pF;
                }else if(mR==0){
                    ed=0;
                    rp=pR;
                }else{
                    ed=mF;
                    rp=pF;
                    if(mR<mF){
                        ed=mR;
                        rp=pR;
                    }
                }

                double acc=(double)(j-ed)/j; //acc:accuracy
                if(acc>=0.90){
                    int sc=j-ed; //sc:score
                    if(f[i]+sc>f[i+j]){
                        f[i+j]=f[i]+sc;
                        g[i+j]=i;
                        grs[i+j]=rp;
                    }
                }
            }
        }
    }

    auto et=high_resolution_clock::now(); //et:end time
    duration<double>df=et-st; //df:difference time
    cout<<"这段代码耗时:"<<df.count()<<"秒\n";
    
    vector<Res>res; //res:final results
    int ce=lQ; //ce:curr_end
    
    while(ce>0){
        int cs=g[ce]; //cs:curr_start
        
        if(cs==ce-1&&grs[ce]==-1){
            ce=cs;
            continue;
        }

        while(cs>0){
            int ps=g[cs]; //ps:prev_start
            if(ps==cs-1&&grs[cs]==-1){
                break;
            }
            
            if(getPM(ps+1,ce)!=-1){
                cs=ps;
            }else{
                break;
            }
        }
        
        int rsp; //rsp:r_start_pos
        if(cs==g[ce]){
            rsp=grs[ce];
        }else{
            rsp=getPM(cs+1,ce);
        }
        
        int qsp=cs; //qsp:q_start_pos
        int qep=ce-1; //qep:q_end_pos
        int rep=rsp+(ce-cs)-1; //rep:r_end_pos
        
        res.push_back({qsp,qep,rsp,rep});
        ce=cs;
    }
    
    reverse(res.begin(),res.end());
    cout<<"[";
    for(size_t i=0;i<res.size();++i){
        cout<<"("<<res[i].qs<<","<<res[i].qe<<","<<res[i].rs<<","<<res[i].re<<")";
        if(i!=res.size()-1){
            cout<<",";
        }
    }
    cout<<"]\n";

    return 0;
}