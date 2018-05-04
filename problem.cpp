#include "problem.h"

Problem::Problem(QWidget *parent)
{

}

Problem::~Problem(){

}

void Problem::Clear(){
    PM.clear();
    PL.clear();
    PIN.clear();
    INDEX.clear();
    TYPE.clear();
    nodes.clear();
    UIN.clear();
    G.clear();
    H.clear();
    A.clear();
    singular=false;
    N=0;
    InnerN=0;
    for(int i=0;i<btables.size();i++){
        btables[i]->setRowCount(0);
        btables[i]->setEditTriggers(QAbstractItemView::AllEditTriggers);
        btables[i]->blockSignals(false);
    }
    innertable->setRowCount(0);
    innertable->setEditTriggers(QAbstractItemView::AllEditTriggers);
    innertable->blockSignals(false);
}

void Problem::AddBoundaryTable(QTableWidget *table){
    btables.push_back(table);
    btables[btables.size()-1]->setColumnCount(6);
    btables[btables.size()-1]->setHorizontalHeaderItem(0,new QTableWidgetItem("X"));
    btables[btables.size()-1]->setHorizontalHeaderItem(1,new QTableWidgetItem("Y"));
    btables[btables.size()-1]->setHorizontalHeaderItem(2,new QTableWidgetItem("Known"));
    btables[btables.size()-1]->setHorizontalHeaderItem(3,new QTableWidgetItem("Type"));
    btables[btables.size()-1]->setHorizontalHeaderItem(4,new QTableWidgetItem("U"));
    btables[btables.size()-1]->setHorizontalHeaderItem(5,new QTableWidgetItem("Q"));
}

void Problem::AddInnerTable(QTableWidget *table){
    innertable=table;
    innertable->setColumnCount(5);
    innertable->setHorizontalHeaderItem(0,new QTableWidgetItem("X"));
    innertable->setHorizontalHeaderItem(1,new QTableWidgetItem("Y"));
    innertable->setHorizontalHeaderItem(2,new QTableWidgetItem("U"));
    innertable->setHorizontalHeaderItem(3,new QTableWidgetItem("Ux"));
    innertable->setHorizontalHeaderItem(4,new QTableWidgetItem("Uy"));

}

void Problem::DeleteBoundaryPoint(int index){
    PL.erase(PL.begin()+index);
    INDEX.erase(INDEX.begin()+index);
    TYPE.erase(INDEX.begin()+index);
    nodes.erase(nodes.begin()+index);
    N--;
    btables[0]->removeRow(index);
    emit UpdateBoundaryPoints(PL);
}

void Problem::DeleteInnerPoint(int index){
    PIN.erase(PIN.begin()+index);
    InnerN--;
    innertable->removeRow(index);
    emit UpdateInnerPoints(PIN);
}

void Problem::AddBoundaryPoint(float x,float y){
    N++;
    PL.push_back(Point(x,y));
    INDEX.push_back(0);
    TYPE.push_back(0);
    nodes.push_back(node());
    btables[0]->blockSignals(true);
    btables[0]->setRowCount(N);
    btables[0]->setItem(N-1,0,new QTableWidgetItem(QString::number(PL[N-1].x)));
    btables[0]->setItem(N-1,1,new QTableWidgetItem(QString::number(PL[N-1].y)));
    QTableWidgetItem *item= new QTableWidgetItem();
    item->setCheckState(Qt::Unchecked);
    btables[0]->setItem(N-1,2,item);
    btables[0]->setItem(N-1,3,new QTableWidgetItem(QString::number(0)));
    btables[0]->setItem(N-1,4,new QTableWidgetItem(QString::number(0)));
    item= new QTableWidgetItem();
    item->setFlags(item->flags()^Qt::ItemIsEditable);
    btables[0]->setItem(N-1,5,item);
    btables[0]->blockSignals(false);
    emit UpdateBoundaryPoints(PL);

}
void Problem::AddInnerPoint(float x,float y){
    InnerN++;
    PIN.push_back(Point(x,y));
    innertable->blockSignals(true);
    innertable->setRowCount(InnerN);
    innertable->setItem(InnerN-1,0,new QTableWidgetItem(QString::number(PIN[InnerN-1].x)));
    innertable->setItem(InnerN-1,1,new QTableWidgetItem(QString::number(PIN[InnerN-1].y)));
    QTableWidgetItem *item= new QTableWidgetItem();
    item->setFlags(item->flags()^Qt::ItemIsEditable);
    innertable->setItem(InnerN-1,2,item);
    QTableWidgetItem *itemx= new QTableWidgetItem();
    itemx->setFlags(itemx->flags()^Qt::ItemIsEditable);
    innertable->setItem(InnerN-1,3,itemx);
    QTableWidgetItem *itey= new QTableWidgetItem();
    itey->setFlags(itey->flags()^Qt::ItemIsEditable);
    innertable->setItem(InnerN-1,4,itey);
    innertable->blockSignals(false);
    emit UpdateInnerPoints(PIN);
}


void Problem::UpdateBoundaryPoint(int index, float x,float y,int Known,int type,float u){
    PL[index].x=x;
    PL[index].y=y;
    INDEX[index]=Known;
    TYPE[index]=type;
    nodes[index].U=u;
    emit UpdateBoundaryPoints(PL);
}

void Problem::UpdateInnerPoint(int index,float x, float y){
    PIN[index].x=x;
    PIN[index].y=y;
    emit UpdateInnerPoints(PIN);
}



void Problem::SaveToFile(QString filename){
    if(PL.size()>0){
        FILE* f;
        f=fopen(filename.toStdString().c_str(),"w");
        fprintf(f,"%d\n",N);
        for(int i=0;i<N;i++){
            fprintf(f,"%f\n%f\n",PL[i].x,PL[i].y);
        }
        for(int i=0;i<N;i++){
            fprintf(f,"%d\n%d\n%f\n",INDEX[i],TYPE[i],nodes[i].U);
        }
        fprintf(f,"%d\n",InnerN);
        for(int i=0;i<InnerN;i++){
            fprintf(f,"%f\n%f\n",PIN[i].x,PIN[i].y);
        }
        fclose(f);
    }
}

void Problem::LoadFromFile(QString filename){
    FILE* f;
    f=fopen(filename.toStdString().c_str(),"r");
    fscanf(f,"%d",&N);
    PL=std::vector<Point>(N);
    for(int i=0;i<N;i++){
        fscanf(f,"%f%f",&PL[i].x,&PL[i].y);
    }
    for(int i=0;i<N;i++){
        int curindex,curType;
        node currentNode=node();
        float currentValue;
        fscanf(f,"%d%d%f",&curindex,&curType,&currentValue);
        INDEX.push_back(curindex);
        TYPE.push_back(curType);
        if(curindex==0){
            if(curType==0){
                currentNode.U=currentValue;
            }
        }
        else{
            currentNode.Q=currentValue;
        }
        nodes.push_back(currentNode);
    }
    fscanf(f,"%d",&InnerN);
    for(int i=0;i<InnerN;i++){
        float x,y;
        fscanf(f,"%f%f",&x,&y);
        PIN.push_back(Point(x,y));
    }
    innertable->blockSignals(true);
    innertable->setRowCount(InnerN);

    for(int i=0;i<InnerN;i++){
        innertable->setItem(i,0,new QTableWidgetItem(QString::number(PIN[i].x)));
        innertable->setItem(i,1,new QTableWidgetItem(QString::number(PIN[i].y)));
        QTableWidgetItem *item= new QTableWidgetItem();
        item->setFlags(item->flags()^Qt::ItemIsEditable);
        innertable->setItem(i,2,item);
    }
    innertable->blockSignals(false);
    btables[0]->blockSignals(true);
    btables[0]->setRowCount(N);

    for(int i=0;i<N;i++){
        btables[0]->setItem(i,0,new QTableWidgetItem(QString::number(PL[i].x)));
        btables[0]->setItem(i,1,new QTableWidgetItem(QString::number(PL[i].y)));
        QTableWidgetItem * item= new QTableWidgetItem();
        if(INDEX[i]==0){
            item->setCheckState(Qt::Checked);
        }
        else{
            item->setCheckState(Qt::Unchecked);
        }
        btables[0]->setItem(i,2,item);
        btables[0]->setItem(i,3,new QTableWidgetItem(QString::number(TYPE[i])));
        btables[0]->setItem(i,4,new QTableWidgetItem(QString::number(nodes[i].U)));
        item= new QTableWidgetItem();
        item->setFlags(item->flags()^Qt::ItemIsEditable);
        btables[0]->setItem(i,5,item);
    }
    btables[0]->blockSignals(false);
    emit UpdateBoundaryPoints(PL);
    emit UpdateInnerPoints(PIN);
}

void Problem::Solve(){
    initialized=false;
    if(N!=0){
        gmatr();
        hmatr();
        int minimumIterations=0;
        do{
        //TODO apo nodes
            fmatr();
            abmatr();
            solveq();
            reorder();// update nodes
            if(singular){
                emit Error("System is Singular "+QString::number(minimumIterations));
            }
            minimumIterations++;
        }while(iterationCheck()||minimumIterations<5);
        uinter();
        deriv();
        UpdateResults();

    }
    else{
        emit Error("Insufficient Data to Solve");
    }
}

void Problem::UpdateResults(){
    btables[0]->blockSignals(true);
    innertable->blockSignals(true);
    for(int i=0;i<N;i++){
        QTableWidgetItem *item= new QTableWidgetItem(QString::number(nodes[i].U));
        btables[0]->setItem(i,4,item);
        item= new QTableWidgetItem(QString::number(nodes[i].Q));
        btables[0]->setItem(i,5,item);
    }

    for(int i=0;i<InnerN;i++){
        QTableWidgetItem *item= new QTableWidgetItem(QString::number(UIN[i]));
        innertable->setItem(i,2,item);
        QTableWidgetItem *itemx= new QTableWidgetItem(QString::number(UPIN[i].x));
        innertable->setItem(i,3,itemx);
        QTableWidgetItem *itemy= new QTableWidgetItem(QString::number(UPIN[i].y));
        innertable->setItem(i,4,itemy);
    }
    btables[0]->setEditTriggers(QAbstractItemView::NoEditTriggers);
    innertable->setEditTriggers(QAbstractItemView::NoEditTriggers);

}

void Problem::gmatr(){

    /* This function computes the elements of the G matrix */




    /* XM and YM: matrices of nodal cordinates */

    PL.push_back(PL[0]);

    for(int i=0; i<N; i++){
        PM.push_back((PL[i]+PL[i+1])/2.0);
    }

    /* Forming G matrix */
    for(int i=0; i<N ; i++){
        G.push_back(std::vector<float>());
        for(int j=0; j<N; j++){


            if(i != j){
                G[i].push_back(rlintc( PM[i],PL[j],PL[j+1]));
            }
            else{
                G[i].push_back(slintc(PL[j],PL[j+1]));
            }
        }
    }
}

void Problem::hmatr(){
    /* XM and YM: matrices of nodal cordinates */
    /*
    // Δεν χρειάζονται πια γιατί έχουν γίνει στην gmatr και διατηρούνται και εδώ.
    PL.push_back(PL[0]);

    for(int i=0; i<N; i++){
        PM.push_back((PL[i]+PL[i])/2.0);
    }
    //*/

    /* Forming H matrix */

    for(int i=0;i<N;i++){
        H.push_back(std::vector<float>());
        for(int j=0;j<N;j++){
            if(i!=j){
                H[i].push_back(dalpha(PM[i],PL[j],PL[j+1]));
            }
            else{
                //H[i].push_back(-0.5);
                H[i].push_back(0.5);
            }

        }
    }

}


void Problem::fmatr(){
    if(initialized){
        for(int i=0;i<N;i++){
            if(TYPE[i]==1){
                // palia sun paragwgos
                F[i]=anode(nodes[i].U);
                qDebug()<<QString::number(i)+" F : "+QString::number(F[i]);
                DF[i]=anodeDF(nodes[i].U);
                qDebug()<<QString::number(i)+" DF : "+QString::number(DF[i]);

            }
            else if(TYPE[i]==2){
                F[i]=cathode(nodes[i].U);
                qDebug()<<QString::number(i)+" F : "+QString::number(F[i]);
                DF[i]=cathodeDF(nodes[i].U);
                qDebug()<<QString::number(i)+" DF : "+QString::number(DF[i]);
            }
        }
    }
    else{
        initialized=true;
        F=std::vector<float>(N,0);
        DF=std::vector<float>(N,0);
        for(int i=0;i<N;i++){
            if(TYPE[i]==1){
                F[i]=anode(initialPhiAnode); //anode  ia*(e^((Φα-Φ)/βα)-e^(-(Φα-Φ)/αα))
                qDebug()<<QString::number(i)+" F : "+QString::number(F[i]);
                DF[i]=anodeDF(nodes[i].U);
                qDebug()<<QString::number(i)+" DF : "+QString::number(DF[i]);
            }
            else if(TYPE[i]==2){
                F[i]=cathode(initialPhiCathode); // cathode
                qDebug()<<QString::number(i)+" F : "+QString::number(F[i]);
                DF[i]=cathodeDF(nodes[i].U);
                qDebug()<<QString::number(i)+" DF : "+QString::number(DF[i]);
            }
        }
    }
}


float Problem::anode(float U){
    return (-Ia/C)*(exp(-(Phia-U)/Aa)-exp((Phia-U)/Ba));
}

float Problem::cathode(float U){
    return (Ic/C)*(exp((U-Phic)/Ac)-exp((Phic-U)/Bc));
}

float Problem::anodeDF(float U){
    return (-Ia/C)*((1/Aa)*exp((-1/Aa)*(Phia-U))+(1/Ba)*exp((1/Ba)*(Phia-U)));
}

float Problem::cathodeDF(float U){
    return (Ic/C)*((1/Bc)*exp((1/Bc)*(Phic-U))+(-1/Ac)*exp((-1/Ac)*(Phic-U)));
}

bool Problem::iterationCheck(){
    bool flag=false;
    for(int i=0;i<N;i++){
        if(INDEX[i]==0&&TYPE[i]>=1){
            if(X[i]>=e){
                flag=true;
            }
        }
    }
    if(flag){
        qDebug("new iteration");
    }
    return flag;
}



void Problem::abmatr(){

    /* This function modifies the G and H matrices by rearranging theirs columns, in order to form
    the A and B=UNB matrices, which must verify the equation [A]{X} = [B] */

    /* A matrix */
    A=std::vector<std::vector<float>>(N,std::vector<float>(N));// orizw tomegethostou pinaka A = NxN
    for(int j=0;j<N;j++){

        /* INDEX(j)=0 --> flux: unknown boundary value ??????????????????????
           INDEX(j)=1 --> potential: unknown boundary ????????????? grammi 175 swsto????????
           UB: the matrix that contains the known boundary conditions*/
        if(INDEX[j]==0){
            switch(TYPE[j]){
            case 0:
                for(int i=0;i<N;i++)
                    A[i][j]=-G[i][j];
                break;
            case 1:
            case 2:
                for(int i=0;i<N;i++)
                    A[i][j]=-G[i][j]*DF[i]+H[i][j];
                break;
            }

        }
        else if(INDEX[j]==1){
            for(int i=0;i<N;i++)
                A[i][j]=H[i][j];
        }
    }

    /* The elements of B matrix will be stored in UNB */

    for(int i=0;i<N;i++){
        B.push_back(0);

        /* INDEX(j)=0 --> potential: known boundary condition
           INDEX(j)=1 --> flux: known boundary contition */

        for(int j=0;j<N;j++){

            if(INDEX[j]==0){
                switch(TYPE[j]){
                case 0:
                    B[i]-= (H[i][j] * nodes[j].U);
                    break;
                case 1:
                case 2:
                    B[i]+=G[i][j]*F[i]-H[i][j]*nodes[i].U;
                    break;
                }

            }
            else if(INDEX[j]==1){
               B[i]+= (G[i][j] * nodes[j].Q);
            }
        }

    }
}

void Problem::solveq(){
    singular=legs();
}


bool Problem::legs(){
    int imax;
    float factor;
    float amax, atemp;
    X=std::vector<float>(N,0);
    for(int j=0; j<N; j++){
        /*
        for(ii=0;ii<N;ii++){
            for(jj=0;jj<N;jj++){
                printf("%f,",a[ii*N+jj]);
            }
            printf("\n");
        }
        //*/
        amax=0.0;
        for(int i=j;i<N; i++){
            if(fabs(amax)<fabs(A[i][j])){
                amax=A[i][j];
                imax=i;
            }
        }
        if(fabs(amax)==0.0){
            return true ;
        }
        //allazw grammh j me imax an den einai isa
        if(j!=imax){
            for(int k=j;k<N; k++){
                atemp=A[j][k];
                A[j][k]=A[imax][k];
                A[imax][k]=atemp;
            }
            atemp=B[imax];
            B[imax]=B[j];
            B[j]=atemp;
        }


        if (j!=N){
            //apaloifh
            for(int i=j+1; i<N; i++){
                factor=A[i][j]/amax;
                for(int jx=j;jx<N; jx++){
                    A[i][jx]-=factor*A[j][jx];
                }
                B[i]-=factor*B[j];
            }
        }
        else{
            break;
        }
    }


    for(int j=N-1;j>=0;j--){
        atemp=0;
        for(int k=j+1;k<N; k++){
            atemp+=A[j][k]*B[k];
        }
        X[j]=(B[j]-atemp)/A[j][j];
    }

    return false;
}


void Problem::reorder(){

    /* This function places the values of the potential u in matrix UB
      and the values of the flux are stored in UNB */
    for(int i=0;i<N;i++){
        if(INDEX[i]==1){
            /* Swap the values of the potential u and the flux */
            nodes[i].U=X[i];
        }
        else if(INDEX[i]==0){
            switch(TYPE[i]){
            case 0:
                nodes[i].Q=X[i];// Gnwsto U
                break;
            case 1:
                nodes[i].U+=X[i];// Anodos
                nodes[i].Q=anode(nodes[i].U);
                break;
            case 2:
                nodes[i].U+=X[i];// Kathodos
                nodes[i].Q=cathode(nodes[i].U);
                break;
            }
        }
    }

}


void Problem::uinter(){

    /* This function computes the values of potential u at the internal points */

    for(int k=0;k<InnerN;k++){
        UIN.push_back(0);
        for(int j=0;j<N;j++){
            float r1=dalpha(PIN[k],PL[j],PL[j+1]);
            float r2=rlintc(PIN[k],PL[j],PL[j+1]);
            //UIN[k]=UIN[k] + r1*UB[j] - r2*UNB[j];
            UIN[k]=UIN[k] - r1*nodes[j].U + r2*nodes[j].Q;
        }
    }
}

void Problem::deriv(){

    /* This function computes the values of the derivatives(fluxes) Ux & Uy at the internal points XIN, YIN */
    float R, SL,angle;
    Point A, B,EN,GP,HP,D,RP,RNT,U,UN;
    /* Gauss integration points and weights */
    Point C[4];
    float xi[4]={-0.86113631, -0.33998104, 0.33998104, 0.86113631};
    float wg[4]={0.34785485, 0.65214515, 0.65214515, 0.34785485};

    for(int i=0;i<InnerN;i++){
        UPIN.push_back(Point());
        for(int j=0;j<N;j++){
            A=(PL[j+1]-PL[j])/2.0;
            B=(PL[j+1]+PL[j])/2.0;
            SL=A.mag();
            angle=atan2(A.y,A.x)-PI/2.0;
            EN=Point(cos(angle),sin(angle));
            GP=Point();
            HP=Point();
            for(int k=0;k<4;k++){
                C[k]=A*xi[k]+B;
                D=C[k]-PIN[i];
                R=D.mag();
                RP=-(D/R);
                RNT=-(RP*EN);
                U=RP/(2*PI*R);
                UN=-Point(RP.x*RNT.x-RP.y*RNT.y,RP.x*RNT.y+RP.y*RNT.x)/(2.0*PI*R*R);
                GP+=U*wg[k]*SL;
                HP+=UN*wg[k]*SL;
            }
            UPIN[i]=UPIN[i]+HP*nodes[j].U-GP*nodes[j].Q;
        }
    }
}




float Problem::rlintc(Point p1,Point p2,Point p3){

    /* This function computes the off-diagonal elements of the G matrix */

    float result, RA, SL;
    Point A,B;
    std::vector<Point> PC(4);

    /* Gauss integration points and weights */
    float xi[4]={-0.86113631, -0.33998104, 0.33998104, 0.86113631};
    float wg[4]={0.34785485, 0.65214515, 0.65214515, 0.34785485};

    A=(p3-p2)/2.0;
    B=(p3+p2)/2.0;

    result=0;

    for(int i=0;i<4;i++){
        PC[i]=A*xi[i]+B;
        RA=(PC[i]-p1).mag();
        //result=result+ log(RA)*wg[i];
        result=result- log(RA)*wg[i];

    }
    //SL=2.0*A.mag();
    //result=result*SL/(4.0*PI);

    SL=A.mag();
    result=result*SL/(2.0*PI);

    return result;
}

float Problem::slintc(Point p1,Point p2){


    /* This function computes the diagonal elements of G matrix */

    float SL, result;
    Point A;
    A=(p2-p1)/2.0;
    SL=A.mag();
    //result=SL*(log(SL)-1.0)/PI;

    result=-SL*(log(SL)-1.0)/PI;

    return result;
}


float Problem::dalpha(Point p1,Point p2,Point p3){

    /* This function computes the off-diagonal elements of the H matrix */

    Point D1,D2,D2R,CosSin;
    float DL1, DA, result;

    D1=p2-p1;
    D2=p3-p1;
    DL1=D1.mag();
    CosSin=D1/DL1;
    D2R=D2*CosSin;
    DA=atan2(D2R.y, D2R.x);
    //result= DA/(2.0*PI);

    result= -DA/(2.0*PI);

    return result;

}
