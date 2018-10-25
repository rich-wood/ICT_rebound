%%%%%%%%%%%%%%%%
% Code for calculation of rebound effects of ICT for Sweden and EU
% Richard Wood, Feb 2017.
% richard.wood@ntnu.no
% Note, running this code requires data structures from EXIOBASE3,
% available with this package or from the author on request.
% Scenario specfication is done externally, and coded here.


%meta information that was pre-compiled at NTNU for EXIOBASE (this based on
%3.3 constant price version)
load(IOmeta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SWE_CNT_INDEX=25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up product names and indices
%%
ict_pro(1).name={'p40.11'};
ict_pro(2).name={'p22'};
ict_pro(3).name={'p30'};
ict_pro(4).name={'p32'};
ict_pro(5).name={'p64'};
ict_pro(6).name={'p72'};
for i=1:6
    ict_pro(i).indx=(strncmp(meta.labsZ(1:200,3),ict_pro(i).name{1},length(ict_pro(i).name{1})));
end

ict_pro(7).name={'ICT_Tot'};
ict_pro(7).indx=logical(sum([ict_pro(3:6).indx],2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up Scenarios
%%

eff(1).indx=ones(200,1);
eff(1).name={['No reduction']};
for i=3:6
    eff(i-1).indx=1-0.1*ict_pro(i).indx; %Assumption 2
    eff(i-1).name={['Reduction in ',ict_pro(i).name{1}]};
end
eff(6).indx=1-0.1*ict_pro(7).indx; %Assumption 2
eff(6).name={['Reduction in ',ict_pro(7).name{1}]};
eff(7).indx=1-0.1*ict_pro(1).indx; %NOTE THIS IS 10% OF TOTAL ELEC CONSUMPTION; WE DON'T HAVE ELEC CONSUMPTION SOLELY FOR ICT.
eff(7).name={['Reduction in ',ict_pro(1).name{1}]};


%%
reb(1).indx=zeros(200,1); %Assumption 0
reb(1).name={['No rebound']};
reb(2).indx=ones(200,1); %Assumption 1
reb(2).name={['Rebound all products']};
for i=3:6
    reb(i).indx=ict_pro(i).indx; %Assumption 2
    reb(i).name={['Rebound in ',ict_pro(i).name{1}]};
end
reb(7).indx=ict_pro(7).indx; %Assumption 2
reb(7).name={['Rebound in ',ict_pro(7).name{1}]};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations
%%

%%


yr_range=1995:2011
% yr_range=2011

%% load the background data
loaddata=0 %only if on NTNU servers, or with access to compiled files to be loaded, otherwise skip, and load pre-loaded multipliers/demand
if loaddata==1
% here we just extract the total footprint and total expenditure by product
% and region from EXIOBASE data, it is saved in Qy_ts.mat
    for yr=yr_range
        yrcnt=yr-min(yr_range)+1
        load(['\\felles.ansatt.ntnu.no\ivt\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_constant_price\Analysis\data\Q',num2str(yr),'.mat'],'Q10','F')
        load(['\\felles.ansatt.ntnu.no\ivt\indecol\Projects\MRIOs\EXIOBASE3\EXIOBASE_3_3_constant_price\Analysis\data\macro',num2str(yr),'.mat'],'macro')
        for reg=1:2
            % total impacts and expenditures for region:
            if reg==1
                y_reg_time(:,yrcnt,reg)=macro.yall.sec.eu(:,1);
                q_reg_time(:,:,yrcnt,reg)=Q10.cons_ind.sec.eu(:,:,1);
            else
                y_reg_time(:,yrcnt,reg)=macro.yall.sec.cnt(:,SWE_CNT_INDEX);
                q_reg_time(:,:,yrcnt,reg)=Q10.cons_ind.sec.cnt(:,:,SWE_CNT_INDEX);
            end
        end
    end
    save('Qy_ts.mat','y_reg_time','q_reg_time','yr_range')
    clear macro Q10 F
else % otherwise load pre-loaded multipliers/demand:
    load('Qy_ts.mat','y_reg_time','q_reg_time','yr_range')
end


%%
disp(' currently respending by average expenditure')
% for marginal in 2010
% load('\\felles.ansatt.ntnu.no\ivt\indecol\Projects\MRIOs\IncomeGrowth\Country Specific Consumption Forecasting Scenario\Results\NoDistribution_LogLog\Engel_Curves_Regressions.mat')
% respend_vector=Engel_Curves.Shares.SE/sum(Engel_Curves.Shares.SE);



for yr=yr_range
    yrcnt=yr-min(yr_range)+1
    for reg=1:2 % we loop over EU and Sweden for results
        
        y_all=y_reg_time(:,yrcnt,reg); % the total HHFD of the region
        q_all=q_reg_time(:,:,yrcnt,reg); % the total Impact by stressor by product of the region
        
        respend_vector=y_all/sum(y_all); % the % of expenditure for rebound calc
        qm_all=q_all./repmat(eps+y_all',size(q_all,1),1); % the multipliers (impact divided by expenditure)
        
        % expenditure and impact of specific product groups
        % this is just the impact that was
        for i=1:7
            %         tmp_F=F.cons_ind.sec.eu(:,ict_pro(i).indx,1);
            tmp_Q=q_all(:,ict_pro(i).indx,1);
            tmp_y=y_all(ict_pro(i).indx,1);
            
            Y(i,yrcnt,reg)=sum(tmp_y);
            QM(i,yrcnt,:,reg)=sum(tmp_Q,2)./repmat(sum(tmp_y),size(tmp_Q,1),1);
            Qres(i,yrcnt,:,reg)=sum(tmp_Q,2);
            %         FM(i,yrcnt,:,reg)=sum(tmp_F,2)./repmat(sum(tmp_y),size(tmp_F,1),1);
        end
        
 %%%%% REBOUND CALCS %%%%%%%%%
        % Now run scenarios
        % Run the scenarios by each efficiency measuer, and for each
        % rebound assumption
        n=0;
        for r=1:7 % efficiency
            for s=1:7 % rebound
                n=n+1;
                y_decreased=y_all.*eff(r).indx;
                savings(n,yrcnt)=sum(y_all)-sum(y_decreased);
                change_respend=(reb(s).indx.*respend_vector)./sum(eps+(reb(s).indx.*respend_vector));
                y_increase=savings(n,yrcnt)*change_respend;
                y_new=y_decreased+y_increase;
                QresSc(n,yrcnt,:,reg)=sum(qm_all.*repmat(y_new',size(q_all,1),1),2); % sum over all products (dim 2)
                QresScName(n,1)={[eff(r).name{1},'; ',reb(s).name{1}]};
                
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Write output
for reg=1:2
    if reg==1
        filename_excel='EU.xls'
    else
        filename_excel='SE.xls'
        copyfile('EU.xls',filename_excel)
    end
    Excel = actxserver ('Excel.Application');
    Excel.Workbooks.Open([pwd,'\',filename_excel]);
    xlswrite1(filename_excel,indic.name,1,'c3')
    xlswrite1(filename_excel,squeeze(Y(:,:,reg)),'Expenditure','d4')
    xlswrite1(filename_excel,[ict_pro.name]','Expenditure','c4')
    xlswrite1(filename_excel,QresScName,'Expenditure','b13')
    xlswrite1(filename_excel,savings,'Expenditure','d13')
    for indx=1:size(QM,3)
        xlswrite1(filename_excel,squeeze(QM(:,:,indx,reg)),num2str(indx),'d4')
        xlswrite1(filename_excel,squeeze(Qres(:,:,indx,reg)),num2str(indx),'d12')
        xlswrite1(filename_excel,indic.name(indx),num2str(indx),'a1')
        xlswrite1(filename_excel,indic.unit(indx),num2str(indx),'a2')
        xlswrite1(filename_excel,QresScName,num2str(indx),'b22')
        xlswrite1(filename_excel,squeeze(QresSc(:,:,indx,reg)),num2str(indx),'d22')
    end
    
    Excel.ActiveWorkbook.Save;
    Excel.ActiveWorkbook.Close;
    Excel.Quit
    clear Excel
end




