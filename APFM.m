% If you use this, please cite 
% Adaptive parallel filter method for active cancellation of road noise inside vehicles

function [w,e,power,mu_save,pos_save] = APFM(I,J,M,Q,B,G,dotnumber,x,d,w_winner,S,mu_f,beta,threshold,delta)

%I ref number, J sec number, M error number, L control filter length, M sec length,dotnumber data length
%B block length, w_winner pre-trained fixed control filter£¬x ref, d disturbance signal
%S sec, mu step size, delta regularization factor, threshold max step size

d = mat2cell(d,ones(1,M));
Nfft=B*2;
sc_frq=cell(J,M); 
power=ones(I,B);
x_f=zeros(I,Nfft);
w_f=cell(J,I);
for j=1:J
    for i=1:I
        w_f{j,i}=fft(w_winner{i,j},Nfft);
    end
end
e_f=cell(M,1);
for k=1:M
    e_f{k,1}=zeros(Nfft,1);
end
vector_x_reverse=zeros(I,B);
vector_x=zeros(I,Nfft);
e=cell(M,1);
e_old=cell(M,Q+1);
vector_e=cell(M,1);
vector_e_old=cell(M,Q+1);
vector_e_all=cell(M,Q+2);

for k=1:M
    e{k,1}=zeros(1,dotnumber);
    vector_e{k,1}=zeros(B,1);
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
        w_old{j,i}(:,1:Q)=w_winner{i,j}(:,1:Q);
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
x_mu_e_f_vector=cell(J,I);
sum_e_all=zeros(1,Q+2);
counter_i=1;
mu_save=zeros(1,dotnumber);
pos=1;
t_block=0;
pos_save=zeros(1,dotnumber);

for n=1:dotnumber
    for i=1:I
        vector_x_reverse(i,:)=[x(i,n),vector_x_reverse(i,1:end-1)];
        vector_x(i,B+counter_i)=x(i,n);
    end
    
    %anti-noise signal   
    for j=1:J
        y_temp=0;
        y_temp_old=zeros(1,Q+1);
        for i=1:I
             y_temp=y_temp+vector_x_reverse(i,:)*w{j,i};
             for q=1:Q+1
                y_temp_old(1,q)=y_temp_old(1,q)+vector_x_reverse(i,:)*w_old{j,i}(:,q);
             end
        end
        y(j)=y_temp;
        vector_y{j}=[y(j);vector_y{j}(1:end-1)];
        for q=1:Q+1
            y_old(j,q)=y_temp_old(1,q);
            vector_y_old{j,q}=[y_old(j,q);vector_y_old{j,q}(1:end-1)];
        end                
    end
	
    %get error signal
    for k=1:M
        e_temp=0;
        e_temp_old=zeros(1,Q+1);
        for j=1:J
            e_temp=e_temp+vector_y{j}.'*S{j,k};
            for q=1:Q+1
                e_temp_old(1,q)=e_temp_old(1,q)+vector_y_old{j,q}.'*S{j,k};
            end
        end
        e{k,1}(1,n)=d{k,1}(1,n)+e_temp;
        vector_e{k}(counter_i)=e{k,1}(1,n);        
        for q=1:Q+1
            e_old{k,q}(1,n)=d{k,1}(1,n)+e_temp_old(1,q);
            vector_e_old{k,q}(counter_i)=e_old{k,q}(1,n);
        end      
    end
    counter_i=counter_i+1; 
	
	%reach a frame
    if counter_i==B+1
        counter_i=1;
        t_block=t_block+1;
        for i=1:I
            x_f(i,:)=fft(vector_x(i,:),Nfft);  
            vector_x(i,1:B)=vector_x(i,B+1:end);
        end
        for k=1:M            
            e_f{k}=fft([zeros(B,1);vector_e{k}],Nfft);
        end
        for j=1:J
            for i=1:I
                x_mu_e_f_vector{j,i}=zeros(Nfft,1);
            end
        end
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
       
        for i=1:I
            for frq_i=1:2*B
                power(i,frq_i)=beta*power(i,frq_i)+(1-beta)*abs(x_f(i,frq_i))*abs(x_f(i,frq_i));

                    for j=1:J
                        n_mu_matrix=0;
                        for k=1:M
                            n_mu_matrix=n_mu_matrix+mu_new_f(i,k)*mu_matrix{j,k}(frq_i)*e_f{k}(frq_i) / power(i,frq_i);
                        end                        
                        filter_x_mul_e_f=conj(x_f(i,frq_i))*n_mu_matrix;
                        x_mu_e_f_vector{j,i}(frq_i)=filter_x_mul_e_f;
                        x_mu_e_f_vector{j,i}(2*B-frq_i+2)=conj( filter_x_mul_e_f);
                    end
            end
        end
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
                sum_ek=sum_ek+sum(vector_e_all{k,q}.^2);
            end
            sum_e_all(1,q)=sum_ek;
        end
        [~,pos]=min(sum_e_all);
                        
        %Delayless
        for j=1:J
            for i=1:I
                w_all{j,i}(:,Q+1)=w{j,i};
                temp=ifft(x_mu_e_f_vector{j,i},Nfft);
                temp(B+1:end)=zeros(B,1);
                delta_w_f=fft(temp);
                if pos==Q+1 %additional initialisation
                    delta_w_f=(1-beta).*delta_w_f;
                    power_x=zeros(1,I);
                    power_e=zeros(1,M);
                    power=ones(I,B);
                end
                w{j,i}=w_all{j,i}(:,pos);
                w_f{j,i}=fft(w{j,i},Nfft);
                w_f{j,i}=w_f{j,i}+delta_w_f;
                w_L=real(ifft(w_f{j,i}));
                w{j,i}=w_L(1:B);
            end
        end    
    end
    pos_save(n)=pos;%save selected filter index   
end
e = cell2mat(e);

end