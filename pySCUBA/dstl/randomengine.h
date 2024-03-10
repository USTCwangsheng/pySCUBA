/*
 * randomengine.h
 *
 *  Created on: 2016年4月21日
 *      Author: hyliu
 */

#ifndef RANDOMENGINE_H_
#define RANDOMENGINE_H_
#include <memory>
#include <utils/randomgenerator.h>

namespace NSPdstl {

template<typename RNGENGINE = std::mt19937>
struct RandomEngine {
public:

	typedef NSPutils::variate_generator<RNGENGINE&, std::uniform_real_distribution<double>> real_generator_type;
	typedef NSPutils::variate_generator<RNGENGINE&, std::uniform_int_distribution<long>> int_generator_type;
	typedef NSPutils::variate_generator<RNGENGINE&, std::normal_distribution<double>> normal_generator_type;

	static RandomEngine &getinstance() {
		static RandomEngine instance;
		if (!instance.initialized_) {
			instance.init();
		}
		return instance;
	}
	template<typename T>
	void init(const T& seed) {
		init();
		reseed(seed);
	}
	void init() {
		engine_ = std::shared_ptr < RNGENGINE > (new RNGENGINE());
		setrealrng(0.0, 1.0);
		initialized_=true;
	}
	template<typename T>
	void reseed(const T & sd) {
		engine_->seed(sd);
	}
	template<typename T1, typename T2>
	void setrealrng(T1 rmin, T2 rmax) {
		std::uniform_real_distribution<double> pdf1(rmin, rmax);
		realrng_.reset();
		realrng_ = std::shared_ptr < real_generator_type
				> (new real_generator_type(*engine_, pdf1));
	}

	template<typename T1, typename T2>
	void setintrng(T1 imin, T2 imax) {
		std::uniform_int_distribution<long> pdf(imin, imax);
		intrng_.reset();
		intrng_ = std::shared_ptr < int_generator_type
				> (new int_generator_type(*engine_, pdf));
	}
	double randomreal() {
		return (*realrng_)();
	}
	long randomint() {
		return (*intrng_)();
	}
	int_generator_type & intrng() {return *intrng_;}
	template<typename T1, typename T2>
	int_generator_type & intrng(T1 imin,T2 imax) {
		setintrng(imin,imax);
		return *intrng_;}

	real_generator_type & realrng() {return *realrng_;}
	template<typename T1, typename T2>
	real_generator_type & realrng(T1 rmin,T2 rmax) {
		setrealrng(rmin,rmax);
		return *realrng_;}
	normal_generator_type & normalrng() {
		if(!normalrng_) {
			std::normal_distribution<double> nd;
			normalrng_=std::shared_ptr<normal_generator_type>
			(new normal_generator_type(*engine_,nd));
		}
		return *normalrng_;
	}
	double randomnormal(double mean=0.0,double sd=1.0){
		return mean+sd*normalrng()();
	}
	std::shared_ptr<RNGENGINE> engine_;
	std::shared_ptr<real_generator_type> realrng_;
	std::shared_ptr<int_generator_type> intrng_;
	std::shared_ptr<normal_generator_type> normalrng_{nullptr};
	RandomEngine() {
		;
	}
	RandomEngine(RandomEngine const&) = delete;
	void operator=(RandomEngine const&) = delete;
private:
	bool initialized_{false};
};
}
#endif /* RANDOMENGINE_H_ */
