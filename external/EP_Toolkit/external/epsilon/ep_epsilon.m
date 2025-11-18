%  [epsilonhat DFn DFd] = ep_epsilon (y, subj, gp, alg, DFn, DFd)
%
%  Repeated-measures ANOVA assumes that error is purely random.
%  A random factor that causes a measurement in one subject
%  to be a higher (or lower) should have no affect on the next
%  measurement in the same subject. This assumption is called
%  circularity or sphericity. Sphericity, is therefore a measure
%  of equality of variances of the differences between treatment
%  levels.
%
%  The adjusted degrees of freedom is returned to correct the p-
%  value obtained from the treatment F statistic:
%
%  padj = 1-fcdf(F,DFn,DFd)
%
%  y is a vector of the data
%  subj is a vector of subject identifiers
%  gp is the number of treatment groups
%  alg is the algorithm used (see below)
%  DFn is the degrees of freedom in the numerator
%  DFd is the degrees of freedom in the numerator
%
%  y values must be ordered appropriately for each subject
%
%  The algorithms include:
%  1) Greenhouse-Geisser (GG)
%     - for 1-way repeated-measures ANOVA (conservative)
%  2) Huynh-Feldt-Lecoutre correction (HFL)
%     - for 1-way repeated-measures ANOVA (less conservative)
%     - for correcting violations of multisample sphericity in
%       repeated-measures designs with >=2 independent groups
%
%  The default is 'HFL', corresponding to the Lecoutre corrected form of
%  the Huynh-Feldt adjustment [3]. The alternative is 'GG' for the
%  Greenhouse-Geisser adjustment.
%
%  Bibliography:
%  [1] Maxwell and Delaney (2004) Designing Experiments and Analyzing
%       Data: A Model Comparison. Psychology Press. Vol 1 p543
%  [2] http://homepages.gold.ac.uk/aphome/spheric.html
%  [3] Keselman, Algina and Kowalchuk (2001) The analysis of repeated
%       measures designs: A review. British Journal of Mathematical
%       and Statistical Psychology (2001), 54, 1-20
%
%  epsilon v1.0 (21/09/2015)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/

function [epsilonhat DFn DFd] = ep_epsilon (y, subj, gp, alg, DFn, DFd)

  % Check input arguments
  if nargin > 6
    error('Invalid number of input arguments')
  end
  if nargin < 4
    alg = 'HFL';
  end

  % Check output arguments
  if nargout > 3
    error('Invalid number of output arguments')
  end

  % Create data matrix M from the input variables
  gnames = unique(subj);
  b = numel(gnames);
  a = numel(y)/b;
  M = zeros(b,a);
  for j=1:b
    if iscell(gnames)
      M(j,:) = y(subj==gnames{j});
    else
      M(j,:) = y(subj==gnames(j));
    end
  end

  % Calculate Greenhouse-Geisser epsilon hat
  E = cov(M);
  mds = mean(diag(E));
  ms = mean(mean(E));
  msr = mean(E,2);
  N = a^2*(mds-ms)^2;
  D = (a-1)*(sum(sum(E.^2))-2*a*sum(msr.^2)+a^2*ms^2);
  epsilonhat = N/D;

  % Original Huynh-Feldt correction
  % Huynh and Feldt (1976) Estimation of the Box correction for degrees of
  % freedom from sample data in randomized block and split-plot designs.
  % Journal of Educational Statistics, 1, 69-82.
  % Where:
  %   b is the number of subjects
  %   a is the number of treatments
  %   gp is the number of subjects
  %
  % (b*(a-1)*epsilonhat-2)/((a-1)*(b-gp-(a-1)*epsilonhat))
  %
  % Huynh-Feldt-Lecoutre correction for better correction for multisample
  % sphericity assumption violations in approximately balanced designs.
  %
  % Lecoutre (1991) A correction for the epsilon Approximate Test in
  % Repeated Measures Designs With Two or More Independent Groups
  % Journal of Educational Statistics. 16(4): 371-372
  % See also Keselman, Algina and Kowalchuk (2001)

  if ~strcmp(alg,'GG') & ~strcmp(alg,'gg')
    epsilonhat = ((b-gp+1)*(a-1)*epsilonhat-2)/((a-1)*(b-gp-(a-1)*epsilonhat));
  end

  % Apply correction to degrees of freedom
  if nargin > 4
    DFn = epsilonhat*DFn;
  else
    DFn = epsilonhat*(a-1);
  end
  if nargin > 5
    DFd = epsilonhat*DFd;
  else
    DFd = epsilonhat*(b-1)*(a-1);
  end
  if epsilonhat < 0.75
    disp('Serious departure from sphericity (epsilonhat < 0.75)');
  end

  % Calculate minimum possible value of epsilon hat
  lowerbound = 1/(a-1);
