/*
 * stringutil.h
 *
 *  Created on: 2022年05月26日
 *      Author: jwzhu
 */

#ifndef RANDOMGENERATOR_H_
#define RANDOMGENERATOR_H_
#include <random>
namespace NSPutils {

template<class Engine, class Distribution>
class variate_generator
{
public:
    typedef Distribution distribution_type;
    typedef typename Distribution::result_type result_type;
    
    /**
     * A random variate generator is used to join a random number
     * generator together with a random number distribution.
     */
    variate_generator(Engine e, Distribution d)
      : _eng(e), _dist(d) { }

    /** Returns: distribution()(engine()) */
    result_type operator()() { return _dist(engine()); }
    /**
     * Returns: distribution()(engine(), value).
     */
    template<class T>
    result_type operator()(const T& value) { return _dist(engine(), value); }

    /**
     * Returns: A reference to the associated uniform random number generator.
     */
    Engine& engine() { return std::ref(_eng); }
    /**
     * Returns: A reference to the associated uniform random number generator.
     */
    const Engine& engine() const { return std::ref(_eng); }

    /** Returns: A reference to the associated \random_distribution. */
    distribution_type& distribution() { return _dist; }
    /**
     * Returns: A reference to the associated random distribution.
     */
    const distribution_type& distribution() const { return _dist; }

    /**
     * Precondition: distribution().min() is well-formed
     *
     * Returns: distribution().min()
     */
    result_type min() const { return (distribution().min)(); }
    /**
     * Precondition: distribution().max() is well-formed
     *
     * Returns: distribution().max()
     */
    result_type max() const { return (distribution().max)(); }

private:
    Engine _eng;
    distribution_type _dist;
};
}

#endif 
