function [codonTimes] = speed_analysis( Genes, number)
importModelCodons();
codonTimes = zeros(1,length(Genes.Codons));
for k =1:length(Genes.Codons) %for each codon in the list 
    for i=1:length(AAA)
        if strcmp(Genes.Codons(k),AAA(i))
            codonTimes(k) = VarName2(i);
            %disp('in here');
        end
    end
end
%plotting_function(codonTimes,'Translation time vs position on sequence',number,Genes)
end