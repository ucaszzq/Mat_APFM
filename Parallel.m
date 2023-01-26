% If you use this, please cite 
% Adaptive parallel filter method for active cancellation of road noise inside vehicles
% Submitted to Mechanical Systems and Signal Processing

function [w,e,power,mu_save,pos_save] = Parallel(I,J,M,Q,B,G,dotnumber,x,d,w_trained,S,mu_f,beta,threshold,delta)

%I ref number, J sec number, M error number, L control filter length, M sec length,dotnumber data length
%B block length, w_trained pre-trained fixed control filter£¬x ref, d disturbance signal
%S sec, mu step size, delta regularization factor, threshold max step size

d = mat2cell(d,ones(1,M));
Nfft=B*2;
sc_frq=cell(J,M); 
power=ones(I,B);
w_f=cell(J,I);
for j=1:J
    for i=1:I
        w_f{j,i}=fft(w_trained{i,j},Nfft);
    end
end
delta_w=zeros(L,1);
e_f=cell(M,1);
for k=1:M
    e_f{k,1}=zeros(Nfft,1);
end
vector_x_reverse=zeros(I,B);
x_block_old=zeros(I,B);
x_block=zeros(I,B);
e=cell(M,1);
e_old=cell(M,Q+1);
vector_e=cell(M,1);
e_block=cell(M,1);
vector_e_old=cell(M,Q+1);
vector_e_all=cell(M,Q+2);
vector_r=cell(I,J,K);
r_block_old=cell(I,J,K);
r_block=cell(I,J,K);
for i=1:I
    for j=1:J
        for k=1:K
            vector_r{i,j,k}=zeros(G,1);
			r_block_old{i,j,k}=zeros(G,1);
			r_block{i,j,k}=zeros(G,1);
        end
    end
end
for k=1:M
    e{k,1}=zeros(1,dotnumber);
    vector_e{k,1}=zeros(B,1);
	e_block{k,1}=zeros(B,1);
    for q=1:Q+1
        e_old{k,q}=ones(1,dotnumber);
        vector_e_old{k,q}=zeros(B,1);
        vector_e_all{k,q}=zeros(B,1);
    end
    vector_e_all{k,Q+2}=zeros(B,1);
end
power_x=zeros(1,I);
power_e=zeros(1,M);
w=cell(J,I);
w_old=cell(J,I);
w_all=cell(J,I);
for j=1:J
    for i=1:I
        w{j,i}=zeros(B,1);
        w_old{j,i}=zeros(B,Q+1);
        w_all{j,i}=zeros(B,Q+2);
        w_old{j,i}(:,1:Q)=w_trained{i,j}(:,1:Q);
        w_all{j,i}(:,1:Q)=w_old{j,i}(:,1:Q);
    end
end
mu_matrix=cell(J,M);
mu_new_f=mu_f.*ones(I,M);
for j=1:J
    for k=1:M
        filter_temp = fft(S{j,k},Nfft);
        sc_frq{j,k} = filter_temp(1:B);
        for frq_i=1:B
            mu_matrix{j,k}(frq_i)=-1*conj( sc_frq{j,k}(frq_i) );
        end
    end
end
vector_y=cell(J,1);
vector_y_old=cell(J,Q+1);
for j=1:J
    vector_y{j,1}=zeros(G,1);
    for q=1:Q+1
        vector_y_old{j,q}=zeros(G,1);
    end
end
y=zeros(J,1);
y_old=zeros(J,Q+1);
sum_e_all=zeros(1,Q+2);
counter_i=1;
mu_save=zeros(1,dotnumber);
pos=1;
pos_save=zeros(1,dotnumber);

for n=1:dotnumber
    for i=1:I
        vector_x_reverse(i,:)=[x(i,n),vector_x_reverse(i,1:end-1)]; %step 1 
    end
    
    %anti-noise signal   
    for j=1:J
        y_temp=0;
        y_temp_old=zeros(1,Q+1);
        for i=1:I
             y_temp=y_temp+vector_x_reverse(i,:)*w{j,i}; %step 2
             for q=1:Q+1
                y_temp_old(1,q)=y_temp_old(1,q)+vector_x_reverse(i,:)*w_old{j,i}(:,q); %step 8_1
             end
        end
        y(j)=y_temp;
        vector_y{j}=[y(j);vector_y{j}(1:end-1)]; %step 3
        for q=1:Q+1
            y_old(j,q)=y_temp_old(1,q);
            vector_y_old{j,q}=[y_old(j,q);vector_y_old{j,q}(1:end-1)]; %step 8_2
        end                
    end
	
    %error signal
    for k=1:M
        e_temp=0;
        e_temp_old=zeros(1,Q+1);
        for j=1:J
            e_temp=e_temp+vector_y{j}.'*S{j,k}; %step 4_1
            for q=1:Q+1
                e_temp_old(1,q)=e_temp_old(1,q)+vector_y_old{j,q}.'*S{j,k};
            end
        end
        e{k,1}(1,n)=d{k,1}(1,n)+e_temp; %step 4_2
        vector_e{k}(counter_i)=e{k,1}(1,n); 
        for q=1:Q+1
            e_old{k,q}(1,n)=d{k,1}(1,n)+e_temp_old(1,q); %step 9 
            vector_e_old{k,q}(counter_i)=e_old{k,q}(1,n); %step 10_1
        end      
    end
    
    for i=1:I
        for j=1:J
            for k=1:M
                xijk=vector_x_reverse(i,:)*S{j,k}; %step 5 
                vector_r{ii,jj,kk}=[ xijk; vector_r{ii,jj,kk}(1:end-1)];  %step 6 
            end
        end
    end  
    counter_i=counter_i+1; 
	
	%reach a frame
    if counter_i==B+1 %step 12
        counter_i=1;
    %variable step size strategy
        for i=1:I
            power_x(i)=beta*power_x(i)+(1-beta)*sum(vector_x_reverse(i,:).^2);
        end
        for k=1:M            
            power_e(k)=beta*power_e(k)+(1-beta)*sum(vector_e{k}.^2);    
        end
        if n>B
            for i=1:I
                for k=1:M 
                    mu_new_f(i,k) = 1/(sqrt(power_x(i))+power_e(k)+delta);
                end
            end
            mu_new_f(mu_new_f>threshold) = threshold;             
        end
				
		x_block_old=x_block; %step 13_1
		x_block=vector_x_reverse; %step 13_2
		e_block=vector_e; %step 13_3
		r_block_old=r_block; %step 13_4
		r_block=vector_r; %step 13_5
 
        for k=1:M
            for q=1:Q
                vector_e_all{k,q}=vector_e_old{k,q};
            end
            vector_e_all{k,Q+1}=vector_e{k};
            vector_e_all{k,Q+2}=vector_e_old{k,Q+1};
        end
        for q=1:Q+2
            sum_ek=0;
            for k=1:M
                sum_ek=sum_ek+sum(vector_e_all{k,q}.^2); %step 10_2
            end
            sum_e_all(1,q)=sum_ek;
        end
        [~,pos]=min(sum_e_all); %step 14

                        
%VSS strategy & Adaptive algorithms (normalised frequency-domain block FxLMS [1], frequency-domain partitioned block FxLMS [2], and subband adaptive filters [3]) 
%can be used to obtain the update of the control filter 'delta_w'.
        for j=1:J
            for i=1:I
                if pos==Q+2 %if select null (Q+2)th filter, additional initialisation
                    delta_w=(1-beta).*delta_w;
                    power_x=zeros(1,I);
                    power_e=zeros(1,M);
                    power=ones(I,B);
                end
                w{j,i}=w_all{j,i}(:,pos); %step 15
                w{j,i}=w{j,i}+delta_w; %step 16
                w_all{j,i}(:,Q+1)=w{j,i};%step 17
            end
        end    
    end
    pos_save(n)=pos;%save selected filter index   
end
e = cell2mat(e);

end

%Reference
%[1] S. Zhang, Y. Wang, H. Guo, C. Yang, X. Wang, N. Liu, A normalized frequency-domain block filtered-x lms algorithm for active vehicle interior noise control, Mech. Syst. and Signal Process. 120 (2019) 150-165.
%[2] J. Lorente, M. Ferrer, M. De Diego, A. Gonzalez, GPU implementation of multichannel adaptive algorithms for local active noise control, IEEE Trans. Audio Speech Lang. Process. 22 (2014) 1624-1635.
%[3] J. Buck, D. Sachau, Active headrests with selective delayless subband adaptive filters in an aircraft cabin, Mech. Syst. and Signal Process. 148 (2021) 107164.