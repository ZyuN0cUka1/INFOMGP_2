#ifndef CONSTRAINTS_HEADER_FILE
#define CONSTRAINTS_HEADER_FILE

using namespace Eigen;
using namespace std;

typedef enum ConstraintType { DISTANCE, COLLISION } ConstraintType;   //You can expand it for more constraints
typedef enum ConstraintEqualityType { EQUALITY, INEQUALITY } ConstraintEqualityType;

//there is such constraints per two variables that are equal. That is, for every attached vertex there are three such constraints for (x,y,z);
class Constraint {
public:

	int m1, m2;											//Two participating meshes (can be the same)  - auxiliary data for users (constraint class shouldn't use that)
	int v1, v2;											//Two vertices from the respective meshes - auxiliary data for users (constraint class shouldn't use that)
	double invMass1, invMass2;							//inverse masses of two bodies
	double refValue;									//Reference values to use in the constraint, when needed (like distance)
	RowVector3d refVector;								//Reference vector when needed (like vector)
	double CRCoeff;										//extra velocity bias
	ConstraintType constraintType;						//The type of the constraint, and will affect the value and the gradient. This SHOULD NOT change after initialization!
	ConstraintEqualityType constraintEqualityType;		//whether the constraint is an equality or an inequality

	Constraint(const ConstraintType _constraintType, const ConstraintEqualityType _constraintEqualityType, const int& _m1, const int& _v1, const int& _m2, const int& _v2, const double& _invMass1, const double& _invMass2, const RowVector3d& _refVector, const double& _refValue, const double& _CRCoeff) :constraintType(_constraintType), constraintEqualityType(_constraintEqualityType), m1(_m1), v1(_v1), m2(_m2), v2(_v2), invMass1(_invMass1), invMass2(_invMass2), refValue(_refValue), CRCoeff(_CRCoeff) {
		refVector = _refVector;
	}

	~Constraint() {}



	//computes the impulse needed for all particles to resolve the velocity constraint, and corrects the velocities accordingly.
	//The velocities are a vector (vCOM1, w1, vCOM2, w2) in both input and output.
	//returns true if constraint was already valid with "currVelocities", and false otherwise (false means there was a correction done)
	//currCOMPositions is a 2x3 matrix, where each row is per one of the sides of the constraints; the rest of the relevant variables are similar, and so should the outputs be resized.
	bool resolveVelocityConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currVertexPositions, const MatrixXd& currCOMVelocities,
		const MatrixXd& currAngularVelocities, const Matrix3d& invInertiaTensor1, const Matrix3d& invInertiaTensor2,
		MatrixXd& correctedCOMVelocities, MatrixXd& correctedAngularVelocities, double tolerance) {

		RowVector3d contactNormal = (currVertexPositions.row(0) - currVertexPositions.row(1)).normalized();
		RowVector3d R1, R2; R1 << currVertexPositions.row(0) - currCOMPositions.row(0); R2 << currVertexPositions.row(1) - currCOMPositions.row(1);
		RowVectorXd constGradient(12); constGradient << contactNormal, R1.cross(contactNormal), -contactNormal, -R2.cross(contactNormal);

		correctedAngularVelocities = currAngularVelocities;
		correctedCOMVelocities = currCOMVelocities;
		RowVectorXd V(12); V << currCOMVelocities.row(0), currAngularVelocities.row(0), currCOMVelocities.row(1), currAngularVelocities.row(1);

		if (constraintEqualityType == EQUALITY && abs(constGradient.dot(V)) <= abs(tolerance))
			return true;
		
		MatrixXd invMassMatrix(MatrixXd::Identity(12, 12));
		invMassMatrix.block<3, 3>(0, 0) = MatrixXd::Identity(3, 3) * invMass1;
		invMassMatrix.block<3, 3>(3, 3) = invInertiaTensor1;
		invMassMatrix.block<3, 3>(6, 6) = MatrixXd::Identity(3, 3) * invMass2;
		invMassMatrix.block<3, 3>(9, 9) = invInertiaTensor2;
		double lambda = -(1 + CRCoeff) * constGradient.dot(V) / (constGradient * invMassMatrix * constGradient.transpose()).sum();
		RowVectorXd correctVector = lambda * invMassMatrix * constGradient.transpose();

		correctedAngularVelocities.block<1, 3>(0, 0) << correctVector.segment<3>(3) + currAngularVelocities.row(0);
		correctedAngularVelocities.block<1, 3>(1, 0) << correctVector.segment<3>(9) + currAngularVelocities.row(1);
		correctedCOMVelocities.block<1, 3>(0, 0) << correctVector.segment<3>(0) + currCOMVelocities.row(0);
		correctedCOMVelocities.block<1, 3>(1, 0) << correctVector.segment<3>(6) + currCOMVelocities.row(1);
		if (constraintEqualityType == INEQUALITY)
			return true;
		return false;
	}

	//projects the position unto the constraint
	//returns true if constraint was already valid with "currPositions"
	bool resolvePositionConstraint(const MatrixXd& currCOMPositions, const MatrixXd& currConstPositions, MatrixXd& correctedCOMPositions, double tolerance) {

		double depth = (currConstPositions.row(0) - currConstPositions.row(1)).norm() - refValue;
		correctedCOMPositions = currCOMPositions;

		if (constraintEqualityType == EQUALITY && abs(depth) <= tolerance)
			return true;

		RowVector3d contactNormal = (currConstPositions.row(0) - currConstPositions.row(1)).normalized();
		RowVectorXd constGradient(6); constGradient << contactNormal, -contactNormal;
		MatrixXd invMassMatrix = MatrixXd::Zero(6, 6);
		invMassMatrix.block<3, 3>(0, 0) = MatrixXd::Identity(3, 3) * invMass1;
		invMassMatrix.block<3, 3>(3, 3) = MatrixXd::Identity(3, 3) * invMass2;

		double lambda = -depth / (constGradient * invMassMatrix * constGradient.transpose()).sum();
		RowVectorXd correctVector = invMassMatrix * constGradient.transpose() * lambda;

		correctedCOMPositions.block<1, 3>(0, 0) << currCOMPositions.row(0) + correctVector.segment<3>(0);
		correctedCOMPositions.block<1, 3>(1, 0) << currCOMPositions.row(1) + correctVector.segment<3>(3);

		if (constraintEqualityType == INEQUALITY)
			return true;
		return false;
	}
};



#endif /* constraints_h */
