%% This code is created for visualizing the topology of 2D stress tensor field
%% Author: Junpeng Wang (junpeng.wang@tum.de)
%% Date: 2021-09-14
%% Refer to paper "Stress Topology Analysis for Porous Infill Optimization" by 
%%	Wang, J., Wu, J. and Westermann, R. (2021, arXiv:2108.09675) for detailed explanations.

function TensorTopologyAnalysis2D(stressfileName)
	%%1. Read Data
	LoadStressField(stressfileName);
	
	%%2. Setup
	Setup();

	%%3. Topology Analysis or PSLs Generation
	figure(1); 
	TensorFieldTopologyAnalysis();
	ShowStressTensorTopology();	
end

function Setup()
	global eleSize_;
	global tracingStepWidth_;
	global relaxedFactor_;
	global degeneratePoints_;
	global majorPSLpool_;
	global minorPSLpool_;
	
	relaxedFactor_ = 1.0;
	degeneratePoints_ = DegeneratePointStruct();
	majorPSLpool_ = PrincipalStressLineStruct();
	minorPSLpool_ = PrincipalStressLineStruct();	
	tracingStepWidth_ = eleSize_;
end

function LoadStressField(fileName)
	global boundingBox_;
	global numNodes_;
	global nodeCoords_;
	global numEles_;
	global eNodMat_;
	global eleState_;
	global nodState_;
	global elementsOnBoundary_;
	global nodesOnBoundary_;
	global cartesianStressField_;
	global loadingCond_; 
	global fixingCond_;
	global eleCentroidList_;

	global nelx_;
	global nely_;
	global carNodMapForward_;
	global voxelizedVolume_;	
	
	%%read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	idx = 1;
	while idx
		idx = idx + 1;
		tmp = fscanf(fid, '%s', 1);
		if strcmp(tmp, 'Resolution:'), idx=0; break; end
		if idx>100, error('Wrong Input!'); end
	end
	tmp = fscanf(fid, '%d %d', [1 2]);
	nelx_ = tmp(1); nely_ = tmp(2); 
	tmp = fscanf(fid, '%s', 1);
	lBound = fscanf(fid, '%f %f', [1 2]);
	tmp = fscanf(fid, '%s', 1);
	uBound = fscanf(fid, '%f %f', [1 2]);
	boundingBox_ = [lBound; uBound];
	tmp = fscanf(fid, '%s', 1); 
	numValidEles = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s', 1);
	validElements = fscanf(fid, '%d', [1, numValidEles])';
	validElements = validElements + 1;		

	%%read cartesian stress field
	tmp = fscanf(fid, '%s %s %s %s', 4);
	numStressFields = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s %s', 2); numLoadedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d %f %f', [3, numLoadedNodes]); 
	tmp(1,:) = tmp(1,:)+1; 
	loadingCond_ = tmp';
	tmp = fscanf(fid, '%s %s', 2); numFixedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d', [1, numFixedNodes]); 
	fixingCond_ = tmp'+1;
	tmp = fscanf(fid, '%s %s', 2); numValidNods = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%f %f %f', [3, numValidNods]);
	cartesianStressField_ = tmp';
	fclose(fid);

	%%recover cartesian mesh
	voxelizedVolume_ = zeros(nelx_*nely_,1);
	voxelizedVolume_(validElements) = 1;
	voxelizedVolume_ = reshape(voxelizedVolume_, nely_, nelx_);
	RecoverCartesianMesh();
	numNod2ElesVec = zeros(numNodes_,1);
	for ii=1:numEles_
		iNodes = eNodMat_(ii,:);
		numNod2ElesVec(iNodes,1) = numNod2ElesVec(iNodes,1)+1;
	end
	nodesOnBoundary_ = find(numNod2ElesVec<4);	
	nodState_ = zeros(numNodes_,1); nodState_(nodesOnBoundary_) = 1;
	eleState_ = 4*ones(1, numEles_);	
	%%%%
	allNodes = zeros(numNodes_,1);
	allNodes(nodesOnBoundary_) = 1;	
	tmp = zeros(numEles_,1);
	for ii=1:4
		tmp = tmp + allNodes(eNodMat_(:,ii));
	end
	elementsOnBoundary_ = find(tmp>0);				
	%%%%
	%% element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2)]/4;
	
	loadingCond_(:,1) = carNodMapForward_(loadingCond_(:,1));
	fixingCond_ = carNodMapForward_(fixingCond_);
end

function RecoverCartesianMesh()	
	global nelx_; global nely_;
	global voxelizedVolume_;
	global boundingBox_;
	global numEles_; global numNodes_; global eleSize_;
	global nodeCoords_; global eNodMat_; 
	global carEleMapBack_; global carEleMapForward_;
	global carNodMapBack_; global carNodMapForward_;
	global meshState_; 
	%    __ x
	%   / 
	%  -y         
	%		4--------3
	%	    |		 |		
	%		|		 |
	%		1--------2
	%	rectangular element		
	eleSize_ = min((boundingBox_(2,:)-boundingBox_(1,:))./[nelx_ nely_]);
	carEleMapBack_ = find(1==voxelizedVolume_);	
	numEles_ = length(carEleMapBack_);
	meshState_ = zeros(nelx_*nely_,1);	
	meshState_(carEleMapBack_) = 1;	
	carEleMapForward_ = zeros(nelx_*nely_,1);	
	carEleMapForward_(carEleMapBack_) = (1:numEles_)';
	nodenrs = reshape(1:(nelx_+1)*(nely_+1), 1+nely_, 1+nelx_);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1)+1, nelx_*nely_, 1);
	eNodMat_ = repmat(eNodVec(carEleMapBack_),1,4);
	tmp = [0 nely_+[1 0] -1];
	for ii=1:4
		eNodMat_(:,ii) = eNodMat_(:,ii) + repmat(tmp(ii), numEles_,1);
	end
	carNodMapBack_ = unique(eNodMat_);
	numNodes_ = length(carNodMapBack_);
	carNodMapForward_ = zeros((nelx_+1)*(nely_+1),1);
	carNodMapForward_(carNodMapBack_) = (1:numNodes_)';	
	for ii=1:4
		eNodMat_(:,ii) = carNodMapForward_(eNodMat_(:,ii));
	end
	nodeCoords_ = zeros((nelx_+1)*(nely_+1),2);

	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nelx_:boundingBox_(2,1);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/nely_:boundingBox_(1,2);		
	nodeCoords_(:,1) = reshape(repmat(xSeed, nely_+1, 1), (nelx_+1)*(nely_+1), 1);
	nodeCoords_(:,2) = repmat(ySeed, 1, nelx_+1)';	
	
	nodeCoords_ = nodeCoords_(carNodMapBack_,:);
end

function [detJ, invJ] = CalcJacobi()
	%% only for 1st-order quad element
	global eNodMat_;
	global nodeCoords_;
	global deShapeFuncs_;
	nEND = 2;
	nEGIP = 4;
	detJ = zeros(nEGIP,1);
	invJ = sparse(nEND*nEGIP,nEND*nEGIP);	
	
	probeEleNods = nodeCoords_(eNodMat_(1,:)',:);
	for kk=1:nEGIP
		Jac = deShapeFuncs_(nEND*(kk-1)+1:nEND*kk,:)*probeEleNods;
		detJ(kk) = det(Jac);
		invJ(nEND*(kk-1)+1:nEND*kk, nEND*(kk-1)+1:nEND*kk) = inv(Jac);		
	end	
end

function N = ShapeFunction(s, t)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes
	s = s(:);
	t = t(:);
	N = zeros(size(s,1), 4);
	N(:,1) = 0.25*(1-s).*(1-t);
	N(:,2) = 0.25*(1+s).*(1-t);
	N(:,3) = 0.25*(1+s).*(1+t);
	N(:,4) = 0.25*(1-s).*(1+t);	
end

function dN = DeShapeFunction(s, t)	
	s = s(:);
	t = t(:);
	dN1ds = -(1-t); dN2ds = 1-t; 	dN3ds = 1+t; dN4ds = -(1+t);
	dN1dt = -(1-s); dN2dt = -(1+s); dN3dt = 1+s; dN4dt = 1-s;
	
	dN = zeros(2*length(s), 4);
	dN(1:2:end,:) = 0.25*[dN1ds dN2ds dN3ds dN4ds];
	dN(2:2:end,:) = 0.25*[dN1dt dN2dt dN3dt dN4dt];	
end

function d2Shape = De2ShapeFunction(s, t)	
	s = s(:);
	t = t(:);
	numCoord = length(s);
	dN1dss = 0; dN2dss = 0; dN3dss = 0; dN4dss = 0;
	dN1dtt = 0; dN2dtt = 0; dN3dtt = 0; dN4dtt = 0;
	dN1dst = 0.25; dN2dst = -0.25; dN3dst = 0.25; dN4dst = -0.25;	
	
	d2Shape = repmat(s,3,4);
	d2Shape(1:3:end,:) = repmat([dN1dss	dN2dss	dN3dss	dN4dss], numCoord, 1);
	d2Shape(2:3:end,:) = repmat([dN1dtt	dN2dtt	dN3dtt	dN4dtt], numCoord, 1);
	d2Shape(3:3:end,:) = repmat([dN1dst	dN2dst	dN3dst	dN4dst], numCoord, 1);
end

function gaussIPs = GaussianIntegral();
	sqrt33 = sqrt(3)/3; 
	gaussIPs = [-1 1 1 -1; -1 -1 1 1]' * sqrt33;
end

function TensorFieldTopologyAnalysis()
	global potentialDegenerateElements_;
	global degeneratePoints_;
	global majorPSLpool_;
	global minorPSLpool_;
	
	%%0. Prepare
	PrepareTA();
	
	%%1. identify potential 'degenerate' elements
	ExtractPotentialDegenerateElements();

	%%2. identify degenerate points
	IdentifyingDegeneratePoints(potentialDegenerateElements_);	
	
	%%3. post-processing degenerate points
	degeneratePoints_ = PostProcessDegeneratePoints();
	
	%%4. get the topology skeletons
	ComputeTopologicalSkeletons();	
	
	%%5. write into PSL pools
	if length(degeneratePoints_)
		for ii=1:length(degeneratePoints_)
			if 0==degeneratePoints_(ii).eleIndex, continue; end
			for jj=1:length(degeneratePoints_(ii).majorSkeletons)
				majorPSLpool_(end+1,1) = degeneratePoints_(ii).majorSkeletons(jj);
			end
			for jj=1:length(degeneratePoints_(ii).minorSkeletons)
				minorPSLpool_(end+1,1) = degeneratePoints_(ii).minorSkeletons(jj);
			end
		end
	end	
end

function PrepareTA()
	global detJ_;
	global invJ_;	
	global deShapeFuncs_;
	
	gaussIPs = GaussianIntegral();
	deShapeFuncs_ = DeShapeFunction(gaussIPs(:,1), gaussIPs(:,2));	
	[detJ_, invJ_] = CalcJacobi();	
end

function ExtractPotentialDegenerateElements()
	global numEles_;
	global eNodMat_;	
	global cartesianStressField_;
	global potentialDegenerateElements_;
	potentialDegenerateElements_ = [];
	for ii=1:numEles_
		eleStress = cartesianStressField_(eNodMat_(ii,:)',:);
		opt = DegenrationMeasure(eleStress);
		if 1==opt, potentialDegenerateElements_(end+1,1) = ii; end			
	end
end

function opt = DegenrationMeasure(tar)
	discriminants = DiscriminantConstraintFuncs(tar);
	v1 = discriminants(:,1);
	
	bool1_1 = v1(1)>0 && v1(2)>0 && v1(3)>0 && v1(4)>0;
	bool1_2 = v1(1)<0 && v1(2)<0 && v1(3)<0 && v1(4)<0;
	bool1 = bool1_1 || bool1_2;
	
	v2 = discriminants(:,2);
	bool2_1 = v2(1)>0 && v2(2)>0 && v2(3)>0 && v2(4)>0;
	bool2_2 = v2(1)<0 && v2(2)<0 && v2(3)<0 && v2(4)<0;
	bool2 = bool2_1 || bool2_2;
	
	if bool1 || bool2, opt = 0; else, opt = 1; end	
end

function discriminants = DiscriminantConstraintFuncs(eleStress)
	discriminants = [eleStress(:,1)-eleStress(:,2), eleStress(:,3)];
end

function IdentifyingDegeneratePoints(potentialElements)	
	global consideredDegenerateElements_; 
	global numConsideredDegenerateElements_; 
	global eNodMat_;
	global cartesianStressField_;
	global paraCoordListDegeEles_;
	
	consideredDegenerateElements_ = potentialElements(:);
	numConsideredDegenerateElements_ = length(consideredDegenerateElements_);	
	paraCoordListDegeEles_ = zeros(numConsideredDegenerateElements_,2);
	for ii=1:numConsideredDegenerateElements_
		iEleStress = cartesianStressField_(eNodMat_(potentialElements(ii),:)',:); 
		v1 = DiscriminantConstraintFuncs(iEleStress);
		[paraCoord, ~, ~, ~] = NewtonIteration(v1, zeros(1,size(v1,2)));
		paraCoordListDegeEles_(ii,:) = paraCoord;
	end
end

function [paraCoordinates, res, opt, index] = NewtonIteration(vtxVec, target)
	%% solving a nonlinear system of equasions by Newton-Rhapson's method
	%%	f1(s,t) = tar1
	%%	f2(s,t) = tar2
	opt = 0;
	normTar = norm(target);
	errThreshold = 1.0e-10; RF = 100*errThreshold;	
	s = -0.0; t = -0.0; maxIts = 150;
	index = 0;
	
	for ii=1:maxIts
		index = index+1;
		c0 = ShapeFunction(s, t)';
		dShape = DeShapeFunction(s, t);
		dns = dShape(1,:)';
		dnt = dShape(2,:)';
		d2Shape = De2ShapeFunction(s, t);
		dnss = d2Shape(1,:)';
		dntt = d2Shape(2,:)';
		dnst = d2Shape(3,:)';
		
		q = vtxVec' * c0;
		dqs = vtxVec' * dns;
		dqt = vtxVec' * dnt;

		dfdv1 = [dqs';dqt'];
		b = dfdv1*(q-target');
		if 0==normTar
			res = norm(q-target');
		else
			res = norm(b);
		end
		if res < errThreshold, break; end			
		
		dfdss = vtxVec'*dnss;
		dfdtt = vtxVec'*dntt;
		dfdst = vtxVec'*dnst;
		A11 = dfdss' * (q-target') + norm(dqs)^2;
		A22 = dfdtt' * (q-target') + norm(dqt)^2;
		A12 = dfdst' * (q-target') + dqs'*dqt;
		A21 = A12;
		A = [A11 A12; A21 A22]; x = A\(-b);		
		s = s + x(1); t = t + x(2);	
	end
	if res <= errThreshold && abs(s)<=RF+1 && abs(t)<=RF+1
		opt = 1;
	end
	paraCoordinates = [s t];
end

function extractedDegeneratePoints = PostProcessDegeneratePoints()
	global consideredDegenerateElements_;
	global paraCoordListDegeEles_;
	global numConsideredDegenerateElements_;
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global thresholdDPE_;
	
	extractedDegeneratePoints = DegeneratePointStruct();
	RF = 1.0e-6; %%relaxation factor
	phyCoordList = [];
	index = 0;
	for ii=1:numConsideredDegenerateElements_	
		paraCoord = paraCoordListDegeEles_(ii,:);
		if abs(paraCoord(1))>RF+1 || abs(paraCoord(2))>RF+1, continue; end
		iDegePot = DegeneratePointStruct();
		iDegePot.eleIndex = consideredDegenerateElements_(ii);
		iDegePot.paraCoord = paraCoord;
		tarEleNodeIndices = eNodMat_(iDegePot.eleIndex,:)';
		iNodeCoord = nodeCoords_(tarEleNodeIndices,:);
		shapeFuncs = ShapeFunction(iDegePot.paraCoord(1), iDegePot.paraCoord(2));
		iDegePot.phyCoord = shapeFuncs*iNodeCoord;
		phyCoordList(end+1,:) = iDegePot.phyCoord;	
		iEleStress = cartesianStressField_(tarEleNodeIndices,:);
		iDegePot.cartesianStress = shapeFuncs*iEleStress;
		ps = ComputePrincipalStress(iDegePot.cartesianStress);
		iDegePot.principalStress = ps([4 1]);
		directDegenerancyMetric = abs(iDegePot.principalStress(1)-iDegePot.principalStress(2)) / ...
				abs(iDegePot.principalStress(1)+iDegePot.principalStress(2));	
		if directDegenerancyMetric>thresholdDPE_, continue; end
		iDegePot.directDegenerancyExtentMetric = directDegenerancyMetric;
		index = index + 1;	
		extractedDegeneratePoints(index,1) = iDegePot;				
	end
	if 0 == extractedDegeneratePoints(1).eleIndex, extractedDegeneratePoints(1) = []; end %% There is no degenerate point
end

function val = DegeneratePointStruct()
	vec = struct('ith', 0, 'length', 0, 'vec', [], 'index', []);
	PSL = PrincipalStressLineStruct();
	val = struct(	...
		'eleIndex',							0,	...
		'paraCoord',						[],	...		
		'phyCoord',							[],	...
		'cartesianStress',					[], ...
		'principalStress',					[],	...
		'directDegenerancyExtentMetric', 	[], ...
		'tangentList',						[],	...
		'stress2phy',						[],	...
		'abcd',								[],	...
		'delta',							0,	...
		'majorSkeletons',					PSL,...
		'minorSkeletons',					PSL,...
		'separatrices',						vec	...
	);
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'ith',						0, 	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'principalStressList',		[] ...
	);	
end

function principalStress = ComputePrincipalStress(cartesianStress)
	principalStress = zeros(size(cartesianStress,1), 1+2+1+2);
	iPS = zeros(1, 6);
	for ii=1:size(cartesianStress,1)
		iCartesianStress = cartesianStress(ii,:);
		A = iCartesianStress([1 3; 3 2]);
		[eigenVec, eigenVal] = eig(A);
		iPS([1 4]) = diag(eigenVal);
		iPS([2 3 5 6]) = reshape(eigenVec,1,4);
		principalStress(ii,:) = iPS;
	end		
end

function ComputeTopologicalSkeletons()
	global boundingBox_;
	global nodeCoords_; 
	global eNodMat_;
	global cartesianStressField_;
	global degeneratePoints_;
	global elementsOnBoundary_;
	global invJ_;
	global tracingStepWidth_;
	
	%%1. get the derivatives of cartesian stresses at degenerate points with respect to the cartesian coordinates
	%%  compute rotational invariant
	for ii=1:length(degeneratePoints_)
		s = degeneratePoints_(ii).paraCoord(1); t = degeneratePoints_(ii).paraCoord(2);
		dShape = DeShapeFunction(s, t);
		eleCoords = nodeCoords_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		eleNodeCartesionStresses = cartesianStressField_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		dN2dPhyC = invJ_(1:2,1:2)*dShape;	
		degeneratePoints_(ii).stress2phy = [(dN2dPhyC(1,:)*eleNodeCartesionStresses )' (dN2dPhyC(2,:)*eleNodeCartesionStresses )' ];
		
		a = (degeneratePoints_(ii).stress2phy(1,1) - degeneratePoints_(ii).stress2phy(2,1))/2; 
		b = (degeneratePoints_(ii).stress2phy(1,2) - degeneratePoints_(ii).stress2phy(2,2))/2; 
		c = degeneratePoints_(ii).stress2phy(3,1); 
		d = degeneratePoints_(ii).stress2phy(3,2);
		degeneratePoints_(ii).abcd = [a b c d];
		degeneratePoints_(ii).delta = a*d - b*c; %% Negative -> Trisector; Positive -> Wedge
		
		rawRoots = roots([d (c+2*b) (2*a-d) -c]);
		degeneratePoints_(ii).tangentList = rawRoots(0==imag(rawRoots));		
	end
	
	%%1.1 Exclude Degenerate Points located in the Boundary Elements
	%% "An Experience-based Solution to Ease up the Numerical Instability Caused by the POTENTIAL  Jaggy Boundary of Cartesian Mesh"
	elementsWithDegeneratePoints = [degeneratePoints_.eleIndex];
	[~, boundaryElementIndicesWithDegeneratePoints] = intersect(elementsWithDegeneratePoints, elementsOnBoundary_);
	degeneratePoints_(boundaryElementIndicesWithDegeneratePoints) = [];
	
	%%2. get the topology skeletons
	stopCond = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/tracingStepWidth_);	
	for ii=1:length(degeneratePoints_)
		degeneratePoints_(ii).majorSkeletons = repmat(degeneratePoints_(ii).majorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		degeneratePoints_(ii).minorSkeletons = repmat(degeneratePoints_(ii).minorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		for jj=1:length(degeneratePoints_(ii).tangentList)
			seed = [degeneratePoints_(ii).eleIndex degeneratePoints_(ii).paraCoord];
			iniDir = [1 degeneratePoints_(ii).tangentList(jj) ]; iniDir = iniDir/norm(iniDir);
			[maxPSL, minPSL] = GeneratePrincipalStressLines(seed, iniDir, stopCond, []);
			
			dis0 = degeneratePoints_(ii).phyCoord - maxPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - maxPSL.phyCoordList(end,:); dis1 = norm(dis1);		
			if dis0 < dis1
				maxPSL.eleIndexList = maxPSL.eleIndexList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.phyCoordList = maxPSL.phyCoordList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.principalStressList = maxPSL.principalStressList(maxPSL.midPointPosition:maxPSL.length,:);
				maxPSL.midPointPosition = 1;
				maxPSL.length = length(maxPSL.eleIndexList);
			else
				maxPSL.eleIndexList = flip(maxPSL.eleIndexList(1:maxPSL.midPointPosition,:));
				maxPSL.phyCoordList = flip(maxPSL.phyCoordList(1:maxPSL.midPointPosition,:));
				maxPSL.principalStressList = flip(maxPSL.principalStressList(1:maxPSL.midPointPosition,:)); 
				maxPSL.midPointPosition = 1;
				maxPSL.length = length(maxPSL.eleIndexList);
			end
			degeneratePoints_(ii).majorSkeletons(jj) = maxPSL;
			
			dis0 = degeneratePoints_(ii).phyCoord - minPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - minPSL.phyCoordList(end,:); dis1 = norm(dis1);
			if dis0 < dis1
				minPSL.eleIndexList = minPSL.eleIndexList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.phyCoordList = minPSL.phyCoordList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.principalStressList = minPSL.principalStressList(minPSL.midPointPosition:minPSL.length,:);
				minPSL.midPointPosition = 1; minPSL.length = length(minPSL.eleIndexList);
			else
				minPSL.eleIndexList = flip(minPSL.eleIndexList(1:minPSL.midPointPosition,:));
				minPSL.phyCoordList = flip(minPSL.phyCoordList(1:minPSL.midPointPosition,:));
				minPSL.principalStressList = flip(minPSL.principalStressList(1:minPSL.midPointPosition,:)); 
				minPSL.midPointPosition = 1; minPSL.length = length(minPSL.eleIndexList);
			end		
			degeneratePoints_(ii).minorSkeletons(jj) = minPSL;	
		end
	end
end

function [majorPSL, minorPSL] = GeneratePrincipalStressLines(initialSeed, iniDir, limiSteps, typePSL)
	majorPSL = PrincipalStressLineStruct();
	minorPSL = PrincipalStressLineStruct();
	
	%%1. Spot the Starting Point	
	[eleIndex, phyCoord, principalStress, opt] = PreparingForTracing(initialSeed);
	if ~opt, return; end
	%%2. Compute PSL(s)
	if ~isempty(iniDir)
		psDir = [5 6];
		majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
		psDir = [2 3];
		minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
	else
		switch typePSL
			case 'MAJOR'
				psDir = [5 6];
				majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);			
			case 'MINOR'
				psDir = [2 3];			
				minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);			
			case 'BOTH'
				psDir = [5 6];
				majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);	
				psDir = [2 3];
				minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);					
		end
	end
end

function iPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps)
	global tracingStepWidth_;
	iPSL = PrincipalStressLineStruct();
	PSLphyCoordList = phyCoord;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;
	%% tracing along first direction (v1)
	nextPoint = phyCoord + tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, eleIndex, psDir, limiSteps);
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	%% tracing along second direction (-v1)			
	nextPoint = phyCoord - tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, -iniDir, eleIndex, psDir, limiSteps);
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);				
	end
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	iPSL.midPointPosition = size(phyCoordList,1)+1;	
	%%2.3 finish Tracing the current major PSL			
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.principalStressList = PSLprincipalStressList;				
end

function [eleIndex, phyCoord, principalStress, opt] = PreparingForTracing(initialSeed)
	global nodeCoords_; global eNodMat_;
	global cartesianStressField_;
	eleIndex = 0;
	phyCoord = 0; 
	principalStress = 0;
	if 3==size(initialSeed,2)
		opt = 1;
		formatedSeed = initialSeed;
		eleIndex = formatedSeed(1,1);
	elseif 2==size(initialSeed,2)
		[eleIndex, paraCoordinates, opt] = FindAdjacentElement(initialSeed);		
		if opt
			formatedSeed = [eleIndex paraCoordinates];
		else
			return;
		end
	else
		error('Wrong Input!');
	end
	NIdx = eNodMat_(eleIndex,:)';
	eleNodeCoords = nodeCoords_(NIdx,:);
	eleCartesianStress = cartesianStressField_(NIdx,:);
	paraCoord = formatedSeed(1, 2:3);
	shapeFuncs = ShapeFunction(paraCoord(1), paraCoord(2));	
	phyCoord = shapeFuncs*eleNodeCoords;						
	interpolatedCartesianStress = shapeFuncs*eleCartesianStress;
	principalStress = ComputePrincipalStress(interpolatedCartesianStress);	
end

function [phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, elementIndex, typePSL, limiSteps)
	%% Tracing the PSL by 2-nd order Runge-Kutta Scheme 
	global eNodMat_;
	global cartesianStressField_;
	global tracingStepWidth_; 
	
	phyCoordList = zeros(limiSteps,2);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,6);
	index = 0;	
	
	intgerScheme = 'RK2'; %% 'RK2', 'EULER'
	switch intgerScheme		
		case 'EULER'
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);		
			while 1==bool1
				index = index + 1; if index > limiSteps, index = index-1; break; end
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
				nextDir = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));
					
				if 0 == AngleTerminationCondition(iniDir, nextDir), index = index-1; break; end	
				iniDir = nextDir;
				phyCoordList(index,:) = nextPoint;
				eleIndexList(index,:) = elementIndex;
				principalStressList(index,:) = principalStress;					
				nextPoint = nextPoint + tracingStepWidth_*iniDir;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
			end				
		case 'RK2'
			%%initialize initial k1 and k2
			k1 = iniDir;
			iniPot = nextPoint - k1*tracingStepWidth_;
			midPot = nextPoint - k1*tracingStepWidth_/2;
			
				
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(midPot);
			if bool1
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
				k2 = DirectionSelecting(k1, principalStress(typePSL), -principalStress(typePSL));
				nextPoint = iniPot + tracingStepWidth_*k2;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				while 1==bool1
					index = index + 1; if index > limiSteps, index = index-1; break; end
					%%k1
					cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
					cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
					principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
					k1 = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));	
					if 0 == AngleTerminationCondition(iniDir, k1), index = index-1; break; end
					%%k2
					midPot = nextPoint + k1*tracingStepWidth_/2;
					[elementIndex2, paraCoordinates2, bool1] = FindAdjacentElement(midPot);
					if ~bool1, index = index-1; break; end
					cartesianStress2 = cartesianStressField_(eNodMat_(elementIndex2,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates2(1), paraCoordinates2(2));
					cartesianStressOnGivenPoint2 = shapeFuncs*cartesianStress2;
					principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
					k2 = DirectionSelecting(k1, principalStress2(typePSL), -principalStress2(typePSL));		
					%%store	
					iniDir = k1;
					phyCoordList(index,:) = nextPoint;
					eleIndexList(index,:) = elementIndex;
					principalStressList(index,:) = principalStress;	
					%%next point
					nextPoint = nextPoint + tracingStepWidth_*k2;		
					[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				end		
			end		
	end
	phyCoordList = phyCoordList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);				
end

function val = AngleTerminationCondition(dirct1, dirct2)
	angle = acos((dirct1*dirct2') / (norm(dirct1)*norm(dirct2)));
	if angle > pi/6
		val = 0;
	else
		val = 1;
	end
end

function targetDirection = DirectionSelecting(originalVec, Vec1, Vec2)
	angle1 = acos(originalVec*Vec1');
	angle2 = acos(originalVec*Vec2');
	if angle1 < angle2
		targetDirection = Vec1;
	else
		targetDirection = Vec2;
	end
end

function [nextElementIndex, paraCoordinates, opt] = FindAdjacentElement(physicalCoordinates)
	global nelx_; 
	global nely_; 
	global eleSize_;
	global nodeCoords_; 
	global eNodMat_;
	global carEleMapForward_;
	global boundingBox_;
	Lbound = boundingBox_(1,:);
	nextElementIndex = 0; paraCoordinates = []; opt = 0;
	
	physicalCoordinates = physicalCoordinates - Lbound;
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/eleSize_);
		if eleX<1 || eleX>nelx_, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/eleSize_);
		if eleY<1 || eleY>nely_, return; end
	end				
	
	tarEle = nely_*(eleX-1)+(nely_-eleY+1);
	nextElementIndex = carEleMapForward_(tarEle);
	if nextElementIndex	
		opt = 1;
		relatedNodes = eNodMat_(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-Lbound;
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / eleSize_ - 1;
	end
end

function ShowStressTensorTopology()
	global numEles_;
	global eNodMat_;
	global nodeCoords_;
	global nodesOnBoundary_;
	global degeneratePoints_;
	
	if length(degeneratePoints_) > 0
		%%2. draw degenerate points and topological skeletons
		for ii=1:length(degeneratePoints_)		
			% 	
			for jj=1:length(degeneratePoints_(ii).majorSkeletons)
				if degeneratePoints_(ii).majorSkeletons(jj).length > 2
					plot(degeneratePoints_(ii).majorSkeletons(jj).phyCoordList(:,1), degeneratePoints_(ii).majorSkeletons(jj).phyCoordList(:,2), ...
						'-', 'color', [1 0 0], 'LineWidth', 3); hold on; %%[252 141 98]/255
				end
			end
			for jj=1:length(degeneratePoints_(ii).minorSkeletons)
				if degeneratePoints_(ii).minorSkeletons(jj).length > 2
					plot(degeneratePoints_(ii).minorSkeletons(jj).phyCoordList(:,1), degeneratePoints_(ii).minorSkeletons(jj).phyCoordList(:,2), ...
						'-', 'color', [0 0 1], 'LineWidth', 3); hold on %%[102 194 165]/255
				end
			end
			if degeneratePoints_(ii).delta > 0
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'sk', 'LineWidth', 4, 'MarkerSize', 20); hold on
			else
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'ok', 'LineWidth', 4, 'MarkerSize', 20); hold on
			end	
			
		end
	end

	%%Show silhouette
	edgeIndices = eNodMat_(:, [1 2 2 1  2 3 3 2  3 4 4 3  4 1 1 4])';
	edgeIndices = reshape(edgeIndices(:), 4, 4*numEles_);	
	tmp = zeros(size(nodeCoords_,1),1); tmp(nodesOnBoundary_) = 1;
	tmp = tmp(edgeIndices); tmp = sum(tmp,1);
	boundaryEleEdges = edgeIndices(:,find(4==tmp));
	xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(boundaryEleEdges);
	yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(boundaryEleEdges);		
	cPatchs = zeros(size(yPatchs));
	hd = patch(xPatchs, yPatchs, cPatchs); hold on;
	set(hd, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	
	axis equal; axis tight; axis off;
end