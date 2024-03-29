#include <cmath>
#include <memory>
#include <fstream>
#include <set>
#include <numeric>
#include <random>
#include <iostream>

#include <CandyPretty/CandyPretty.h>

#include <boost/exception/all.hpp>
#include <boost/variant.hpp>

#if 0

/*
                D(t)V(t) = E~(D(T)V(T)|F(T))
 */

/*
        This is the analytic solution to geometric brownian motion with drift
                        
                dS = a S dt + b S dw

 */



#include <ql/pricingengines/blackcalculator.hpp>




struct Differential{
        virtual ~Differential()=default;
        /*
                dx = f(x,dt,dw)
                x(t + dt ) = x(t) + dx(t)
         */
        virtual double Eval(double x, double dt, double std_norm)const=0;
};

struct ProcessContext;

struct ProcessIntegral{
        ProcessIntegral()=default;
        ProcessIntegral(ProcessContext& ctx, double x, std::shared_ptr<Differential> dx);
        void SmallChange(double dt, double std_norm){
                x_ += dx_->Eval( x_, dt, std_norm );
        }
        double Value()const{ return x_; }
private:
        double x_;
        std::shared_ptr<Differential> dx_;
};

struct ProcessContext{
        void Register(ProcessIntegral* ptr){
                procs_.push_back(ptr);
        }
        void Step(double dt){
                auto std_norm = [&](){ return D_(G_); };
                for(auto ptr : procs_){
                        ptr->SmallChange(dt, std_norm());
                }
        }
private:
        #if 0
        std::default_random_engine G_;
        #else
        std::random_device G_;
        #endif
        std::normal_distribution<double> D_{0.0, 1.0};
        std::vector<ProcessIntegral*> procs_;
};

inline ProcessIntegral::ProcessIntegral(ProcessContext& ctx, double x, std::shared_ptr<Differential> dx)
        :x_(x),
        dx_(dx)
{
        ctx.Register(this);
}

struct ProcessView{
        struct Impl{
                virtual ~Impl()=default;
                virtual double Value()const=0;
        };
        struct SptrImpl : Impl{
                explicit SptrImpl(std::shared_ptr<ProcessIntegral const> q_) :q(q_) {}
                virtual double Value()const{ return q->Value(); }
                std::shared_ptr<ProcessIntegral const> q;
        };

        ProcessView()=default;
        ProcessView(std::shared_ptr<ProcessIntegral> p){
                impl_ = std::make_shared<SptrImpl>(p);
        }
        // assumes objects lifetime exists for as long as it'self
        ProcessView(ProcessIntegral const& p){
                std::shared_ptr<ProcessIntegral const> aux(&p, [](auto*){});
                impl_ = std::make_shared<SptrImpl>(aux);
        }
        
        ProcessView& operator=(std::shared_ptr<ProcessIntegral> p){
                impl_ = std::make_shared<SptrImpl>(p);
                return *this;
        }
        ProcessView& operator=(ProcessIntegral const& p){
                std::shared_ptr<ProcessIntegral const> aux(&p, [](auto*){});
                impl_ = std::make_shared<SptrImpl>(aux);
                return *this;
        }
        
        
        double Value()const{ return impl_->Value(); }

        // for printing to csv etc
        // this is how it is set, p.Name()  = "Discount ProcessIntegral()"
        std::string& Name(){ return name_; }
        std::string const& Name()const{ return name_; }
protected:
        std::shared_ptr<Impl> impl_;
        std::string name_{"ProcessIntegral"};
};



struct IdentityDifferential : Differential{
        virtual double Eval(double x, double dt, double std_norm)const override{
                return dt;
        }
};

struct BankAccountDifferential : Differential{
        explicit BankAccountDifferential(ProcessView interest_rate)
                :interest_rate_(interest_rate)
        {}
        virtual double Eval(double x, double dt, double std_norm)const override{
                return x * interest_rate_.Value() * dt;
        }
private:
        ProcessView interest_rate_;
};

struct GeometricBrownianMotionWithDriftDifferential : Differential{
        GeometricBrownianMotionWithDriftDifferential(double S0, double r, double sigma)
                :S0_(S0),
                r_(r),
                sigma_(sigma)
        {}
        virtual double Eval(double x, double dt, double std_norm)const override{
                double a = r_ * dt + sigma_ * std_norm * std::sqrt(dt);
                double b = a * x;
                return b;
        }
private:
        double S0_;
        double r_;
        double sigma_;
};

struct VasicekDifferential : Differential{
        VasicekDifferential(double alpha, double beta, double sigma)
                :alpha_(alpha),
                beta_(beta),
                sigma_(sigma)
        {}
        virtual double Eval(double x, double dt, double std_norm)const override{
                double a = ( alpha_ - beta_ * x ) * dt;
                double b = sigma_ * std_norm * std::sqrt(dt);
                double c = a + b;
                return c;
        }
private:
        double alpha_;
        double beta_;
        double sigma_;
};
struct CoxIngersollRos : Differential{
        CoxIngersollRos(double alpha, double beta, double sigma)
                :alpha_(alpha),
                beta_(beta),
                sigma_(sigma)
        {}
        virtual double Eval(double x, double dt, double std_norm)const override{
                double a = ( alpha_ - beta_ * x ) * dt;
                double b = sigma_ * std::sqrt(x) *  std_norm * std::sqrt(dt);
                double c = a + b;
                return c;
        }
private:
        double alpha_;
        double beta_;
        double sigma_;
};





struct AnaBlack : ProcessView{
        AnaBlack(ProcessView t){
                struct AnaBlackImpl : Impl{
                        AnaBlackImpl(ProcessView t)
                                :t_(t)
                        {}
                        virtual double Value()const override{
                                double r = 0.02;
                                double vol = 0.1;
                                double s0 = 10.0;
                                double k = 1.5 * s0;


                                double T = t_.Value();
                                auto discount = std::exp( -r * T );
                                auto fwd = s0 / discount;
                                auto std_dev = std::sqrt( vol * vol * T);
                                QuantLib::BlackCalculator bc(QuantLib::Option::Call,
                                                             k,
                                                             fwd,
                                                             std_dev,
                                                             discount);
                                return bc.value();
                        }
                private:
                        ProcessView t_;
                };
                impl_ = std::make_shared<AnaBlackImpl>(t);
        }
};

struct DiscountProcess : ProcessView{
        DiscountProcess(ProcessView t, double r){
                struct DPImpl : Impl{
                        DPImpl(ProcessView t_, double r_)
                                :t(t_),
                                r(r_)
                        {}
                        virtual double Value()const override{
                                return std::exp( - t.Value() * r );
                        }
                private:
                        ProcessView t;
                        double r;
                };
                impl_ = std::make_shared<DPImpl>(t,r);
        }
};

struct AverageView : ProcessView{
        struct Final : Impl{
                virtual double Value()const override{
                        size_t n = v_.size();
                        double sigma = 0.0;
                        for(size_t idx=0;idx!=n;++idx){
                                sigma +=  v_[idx].Value();
                        }
                        return sigma / n;
                }
        private:
                friend struct AverageView;
                std::vector<ProcessView> v_;
        };
        AverageView(){
                impl_ = std::make_shared<Final>();
        }
        AverageView& Add(ProcessView view){
                auto casted = dynamic_cast<Final*>(impl_.get());
                casted->v_.push_back(view);
                return *this;
        }
        template<class Iter>
        AverageView(Iter first, Iter last){
                impl_ = std::make_shared<Final>();
                for(;first!=last;++first){
                        Add(*first);
                }
        }
};

struct Option : ProcessView{
        Option(ProcessView process, double strike){
                struct OptionImpl : Impl{
                        OptionImpl(ProcessView process, double strike)
                                :process_(process),
                                strike_(strike)
                        {}
                        virtual double Value()const override{
                                return (std::max)(process_.Value() - strike_, 0.0);
                        }
                private:
                        ProcessView process_;
                        double strike_;
                };
                impl_ = std::make_shared<OptionImpl>(process, strike);
        }
};

struct ProcessViewRenderer{
        ProcessViewRenderer(std::ostream& out, std::vector<ProcessView> const& views)
                :out_{std::shared_ptr<std::ostream>(&out, [](auto*){})}, views_(views)
        {
                EmitHeader_();
        }
        ProcessViewRenderer(std::shared_ptr<std::ostream> out, std::vector<ProcessView> const& views)
                :out_{out}, views_(views)
        {
                EmitHeader_();
        }
        void RenderLine(){
                std::vector<std::string> line;
                for(auto const& view : views_){
                        line.push_back(boost::lexical_cast<std::string>(view.Value()));
                }
                lines_.push_back(std::move(line));
        }
        void Emit(){
                CandyPretty::RenderTablePretty(*out_, lines_, opts_);
        }
private:
        void EmitHeader_(){
                std::vector<std::string> line;
                for(auto const& view : views_){
                        line.push_back(boost::lexical_cast<std::string>(view.Name()));
                }
                lines_.push_back(std::move(line));
        }
        CandyPretty::RenderOptions opts_{CandyPretty::RenderOptions::CsvOptions()};
        std::shared_ptr<std::ostream> out_;
        std::vector<ProcessView> views_;
        std::vector<CandyPretty::LineItem> lines_;
};

void example_0(){
        using namespace CandyPretty;

        double r = 0.02;
        double vol = 0.1;
        double T = 20;
        double s0 = 10.0;

        auto gbm = std::make_shared<GeometricBrownianMotionWithDriftDifferential>(s0, r, vol);

        enum{ SampleSize = 4000 };
        ProcessContext ctx;
        auto t = std::make_shared<ProcessIntegral>(ctx, 0, std::make_shared<IdentityDifferential>() );

        std::vector<std::shared_ptr<ProcessIntegral> > gbm_sample(SampleSize);
        for(size_t idx=0;idx!=SampleSize;++idx){
                gbm_sample[idx] = std::make_shared<ProcessIntegral>(ctx, s0, gbm);
        }

        DiscountProcess disc(t, r);

        std::vector<AverageView> avg_vec;
        auto first = &gbm_sample[0];
        avg_vec.emplace_back( first, first + 10);
        avg_vec.back().Name() = "Avg_{10}";
        avg_vec.emplace_back( first, first + 100);
        avg_vec.back().Name() = "Avg_{100}";
        avg_vec.emplace_back( first, first + 1000);
        avg_vec.back().Name() = "Avg_{1000}";
        avg_vec.emplace_back( first, first + 2000);
        avg_vec.back().Name() = "Avg_{2000}";
        avg_vec.emplace_back( first, first + 4000);
        avg_vec.back().Name() = "Avg_{4000}";

        
        std::vector<ProcessView> views;
        views.push_back(t);
        views.back().Name() = "t";
        views.push_back(disc);
        views.back().Name() = "D(t)";
        for(auto const& _ : avg_vec){
                views.push_back(_);
        }
        enum{ GbmViews = 20 };
        for(size_t idx=0;idx < gbm_sample.size() && idx < GbmViews;++idx){
                views.push_back(gbm_sample[idx]);
                std::stringstream sstr;
                sstr << "W_{" << idx << "}(t)";
                views.back().Name() = sstr.str();
        }

        std::vector<ProcessView>* render_view = &views;
        
        #if 0

        AnaBlack black( t);

        std::vector<ProcessView> call_view;
        call_view.push_back(t);
        call_view.push_back(disc);
        call_view.push_back(black);

        std::vector<ProcessView> call_options(SampleSize);
        std::vector<std::pair<AverageView, size_t> > call_avg;
        call_avg.emplace_back(AverageView(), 10);
        call_avg.emplace_back(AverageView(), 100);
        call_avg.emplace_back(AverageView(), 1000);
        call_avg.emplace_back(AverageView(), 2000);
        call_avg.emplace_back(AverageView(), 4000);
        for(auto& _ : call_avg){
                call_view.push_back(_.first);
        }
        for(size_t idx=0;idx!=SampleSize;++idx){
                auto& p = gbm_sample[idx];
                call_options[idx] = Option(p, k);
                for(auto& _ : call_avg){
                        if( idx < _.second ){
                                _.first.Add(call_options[idx]);
                        }
                }
        }







        double ir_0 = 0.05;
        auto f = 10.0;
        auto ir_diff = std::make_shared<VasicekDifferential>(1/f, 20/f, 0.1);
        std::vector<ProcessIntegral> bank_acct_sample(SampleSize);
        std::vector<ProcessView> ir_view;
        AverageView bank_acct_avg;
        AverageView ir_avg;

        ir_view.push_back(t);
        ir_view.push_back(ir_avg);
        ir_view.push_back(bank_acct_avg);

        for(size_t idx=0;idx!=SampleSize;++idx){
                auto irp_p = std::make_shared<ProcessIntegral>(ir_0, ir_diff);
                ProcessView irp(irp_p);
                procs.push_back(irp_p.get());
                auto bank_acct_diff = std::make_shared<BankAccountDifferential>(irp);
                auto bank_acct_p = std::make_shared<ProcessIntegral>(1.0, bank_acct_diff);
                ProcessView bank_acct(bank_acct_p);
                procs.push_back(bank_acct_p.get());

                if( idx < GbmViews ){
                        ir_view.push_back(irp);
                        ir_view.push_back( bank_acct);
                }
                ir_avg.Add(irp);
                bank_acct_avg.Add(bank_acct);
        }
        #endif
        
        


        size_t N = 1000;
        double dt = T / N;

        std::ofstream of{"RiskNeutralBrownianMotion.csv"};
        if( ! of.is_open() )
                BOOST_THROW_EXCEPTION(std::domain_error("unable to open RiskNeutralBrownianMotion.csv"));
        ProcessViewRenderer renderer{of, *render_view};

        for(size_t idx=0;idx!=N;++idx){
                ctx.Step(dt);
                renderer.RenderLine();
        }
        renderer.Emit();
}


void example_1(){
        using namespace CandyPretty;

        double r = 0.02;
        double vol = 0.1;
        double T = 40;
        double s0 = 10.0;
        double k = 1.5 * s0;

        auto gbm = std::make_shared<GeometricBrownianMotionWithDriftDifferential>(s0, r, vol);

        enum{ SampleSize = 4000 };
        ProcessContext ctx;
        auto t = std::make_shared<ProcessIntegral>(ctx, 0, std::make_shared<IdentityDifferential>() );

        std::vector<std::shared_ptr<ProcessIntegral> > gbm_sample(SampleSize);
        std::vector<ProcessView> call_options(SampleSize);
        for(size_t idx=0;idx!=SampleSize;++idx){
                gbm_sample[idx] = std::make_shared<ProcessIntegral>(ctx, s0, gbm);
                call_options[idx] = Option(gbm_sample[idx], k);
        }

        DiscountProcess disc(t, r);
        
        AnaBlack black( t);

        std::vector<AverageView> avg_vec;
        auto first = &call_options[0];
        avg_vec.emplace_back( first, first + 10);
        avg_vec.back().Name() = "Avg_{10}";
        avg_vec.emplace_back( first, first + 100);
        avg_vec.back().Name() = "Avg_{100}";
        avg_vec.emplace_back( first, first + 1000);
        avg_vec.back().Name() = "Avg_{1000}";
        avg_vec.emplace_back( first, first + 2000);
        avg_vec.back().Name() = "Avg_{2000}";
        avg_vec.emplace_back( first, first + 4000);
        avg_vec.back().Name() = "Avg_{4000}";

        
        std::vector<ProcessView> views;
        views.push_back(t);
        views.back().Name() = "t";
        views.push_back(disc);
        views.back().Name() = "D(t)";
        views.push_back(black);
        views.back().Name() = "BS(.)";
        for(auto const& _ : avg_vec){
                views.push_back(_);
        }
        enum{ GbmViews = 20 };
        for(size_t idx=0;idx < call_options.size() && idx < GbmViews;++idx){
                views.push_back(call_options[idx]);
                std::stringstream sstr;
                sstr << "Call_{" << idx << "}(t)";
                views.back().Name() = sstr.str();
        }

        std::vector<ProcessView>* render_view = &views;
        
        #if 0


        std::vector<ProcessView> call_view;
        call_view.push_back(t);
        call_view.push_back(disc);
        call_view.push_back(black);

        std::vector<ProcessView> call_options(SampleSize);
        std::vector<std::pair<AverageView, size_t> > call_avg;
        call_avg.emplace_back(AverageView(), 10);
        call_avg.emplace_back(AverageView(), 100);
        call_avg.emplace_back(AverageView(), 1000);
        call_avg.emplace_back(AverageView(), 2000);
        call_avg.emplace_back(AverageView(), 4000);
        for(auto& _ : call_avg){
                call_view.push_back(_.first);
        }
        for(size_t idx=0;idx!=SampleSize;++idx){
                auto& p = gbm_sample[idx];
                call_options[idx] = Option(p, k);
                for(auto& _ : call_avg){
                        if( idx < _.second ){
                                _.first.Add(call_options[idx]);
                        }
                }
        }







        double ir_0 = 0.05;
        auto f = 10.0;
        auto ir_diff = std::make_shared<VasicekDifferential>(1/f, 20/f, 0.1);
        std::vector<ProcessIntegral> bank_acct_sample(SampleSize);
        std::vector<ProcessView> ir_view;
        AverageView bank_acct_avg;
        AverageView ir_avg;

        ir_view.push_back(t);
        ir_view.push_back(ir_avg);
        ir_view.push_back(bank_acct_avg);

        for(size_t idx=0;idx!=SampleSize;++idx){
                auto irp_p = std::make_shared<ProcessIntegral>(ir_0, ir_diff);
                ProcessView irp(irp_p);
                procs.push_back(irp_p.get());
                auto bank_acct_diff = std::make_shared<BankAccountDifferential>(irp);
                auto bank_acct_p = std::make_shared<ProcessIntegral>(1.0, bank_acct_diff);
                ProcessView bank_acct(bank_acct_p);
                procs.push_back(bank_acct_p.get());

                if( idx < GbmViews ){
                        ir_view.push_back(irp);
                        ir_view.push_back( bank_acct);
                }
                ir_avg.Add(irp);
                bank_acct_avg.Add(bank_acct);
        }
        #endif
        
        


        size_t N = 1000;
        double dt = T / N;

        std::ofstream of{"CallOption.csv"};
        if( ! of.is_open() )
                BOOST_THROW_EXCEPTION(std::domain_error("unable to open CallOption.csv"));
        ProcessViewRenderer renderer{of, *render_view};

        for(size_t idx=0;idx!=N;++idx){
                ctx.Step(dt);
                renderer.RenderLine();
        }
        renderer.Emit();
}


void example_2(){
        using namespace CandyPretty;

        double r = 0.02;
        double vol = 0.1;
        double T = 40;
        double s0 = 10.0;

        enum{ SampleSize = 4000 };

        ProcessContext ctx;

        auto t = std::make_shared<ProcessIntegral>(ctx, 0, std::make_shared<IdentityDifferential>() );

        std::vector<ProcessView> interest_rate_samples(SampleSize);
        std::vector<ProcessView> bank_account_samples(SampleSize);

        double ir_0 = 0.05;
        auto f = 10.0;
        auto ir_diff = std::make_shared<VasicekDifferential>(1/f, 20/f, 0.1);
        
        for(size_t idx=0;idx!=SampleSize;++idx){
                interest_rate_samples[idx] = std::make_shared<ProcessIntegral>(ctx, ir_0, ir_diff);
                auto bank_acct_diff = std::make_shared<BankAccountDifferential>(interest_rate_samples[idx]);
                bank_account_samples[idx]  = std::make_shared<ProcessIntegral>(ctx, 1.0, bank_acct_diff);
        }

        
        std::vector<AverageView> avg_vec;
        auto first = &bank_account_samples[0];
        avg_vec.emplace_back( first, first + 10);
        avg_vec.back().Name() = "Avg_{10}";
        avg_vec.emplace_back( first, first + 100);
        avg_vec.back().Name() = "Avg_{100}";
        avg_vec.emplace_back( first, first + 1000);
        avg_vec.back().Name() = "Avg_{1000}";
        avg_vec.emplace_back( first, first + 2000);
        avg_vec.back().Name() = "Avg_{2000}";
        avg_vec.emplace_back( first, first + 4000);
        avg_vec.back().Name() = "Avg_{4000}";

        
        std::vector<ProcessView> views;
        views.push_back(t);
        views.back().Name() = "t";
        for(auto const& _ : avg_vec){
                views.push_back(_);
        }
        enum{ GbmViews = 20 };
        for(size_t idx=0;idx < SampleSize && idx < GbmViews;++idx){
                std::stringstream sstr;
                views.push_back(interest_rate_samples[idx]);
                sstr << "R_{" << idx << "}(t)";
                views.back().Name() = sstr.str();
                views.push_back(bank_account_samples[idx]);
                sstr.str("");
                sstr << "B_{" << idx << "}(t)";
                views.back().Name() = sstr.str();
        }

        std::vector<ProcessView>* render_view = &views;
        


        size_t N = 1000;
        double dt = T / N;

        std::ofstream of{"BankAccount.csv"};
        if( ! of.is_open() )
                BOOST_THROW_EXCEPTION(std::domain_error("unable to open BankAccount.csv"));
        ProcessViewRenderer renderer{of, *render_view};

        for(size_t idx=0;idx!=N;++idx){
                ctx.Step(dt);
                renderer.RenderLine();
        }
        renderer.Emit();
}

#endif

struct Omega{};
struct Nul{};
struct Not;
struct Union;
struct Intersection;
struct Interval;

using BorelSet = boost::variant<
        Omega,
        Nul,
        boost::recursive_wrapper<Not>,
        boost::recursive_wrapper<Union>,
        boost::recursive_wrapper<Intersection>,
        boost::recursive_wrapper<Interval>
>;

struct Not{
        BorelSet child;
};
struct Union{
        template<class... Args>
        Union(Args&&... args):children{args...}{}

        void Add(Union const& that){
                for(auto const& _ : that.children)
                        children.push_back(_);
        }

        std::vector<BorelSet> children;
};
struct Intersection{
        template<class... Args>
        Intersection(Args&&... args):children{args...}{}

        void Add(Intersection const& that){
                for(auto const& _ : that.children)
                        children.push_back(_);
        }

        std::vector<BorelSet> children;
};
struct IntervalEndPoint{
        friend std::ostream& operator<<(std::ostream& ostr, IntervalEndPoint const& self){
                ostr << "IntervalEndPoint{is_open = " << self.is_open;
                ostr << ", point = " << self.point << "}";
                return ostr;
        }
        bool operator<(IntervalEndPoint const& that)const{
                if( point != that.point )
                        return point < that.point;
                return is_open < that.is_open;
        }
        bool operator==(IntervalEndPoint const& that)const{
                return point == that.point && is_open == that.is_open;
        }
        bool operator!=(IntervalEndPoint const& that)const{
                return ! operator==(that);
        }
        bool is_open;
        double point;

        IntervalEndPoint Switch()const{ return IntervalEndPoint{ ! is_open, point }; }
};
IntervalEndPoint Open(double x){
        return IntervalEndPoint{true, x};
}
IntervalEndPoint Closed(double x){
        return IntervalEndPoint{false, x};
}

struct IntervalUnion{
        template<class... Args>
        IntervalUnion(Args&&... args):children{args...}{}
        
        void Add(IntervalUnion const& that){
                for(auto const& _ : that.children)
                        children.push_back(_);
        }

        std::vector<Interval> children;

        operator Union()const{
                return AsUnion();
        }
        Union AsUnion()const{
                Union tmp;
                for(auto const& _ : children )
                        tmp.children.push_back(_);
                return tmp;
        }

        #if 0
        bool operator==(IntervalUnion const& that)const{
                return children == that.children;
        }
        bool operator!=(IntervalUnion const& that)const{
                return ! operator==(that);
        }
        #endif
        bool operator<(IntervalUnion const& that)const{
                if( children.size() != that.children.size() ){
                        return  children.size() < that.children.size();
                }
                return children < that.children;
        }
};

struct Interval{
        IntervalEndPoint left;
        IntervalEndPoint right;

        static Interval Closed(double x, double y){
                return Interval{ IntervalEndPoint{ false, x},
                                 IntervalEndPoint{ false, y} };
        }
        static Interval Open(double x, double y){
                return Interval{ IntervalEndPoint{ false, x},
                                 IntervalEndPoint{ false, y} };
        }

        bool operator<(Interval const& that)const{
                if( left != that.left )
                        return left < that.left;
                return right < that.right;
        }


        // homogenous operations are here
        IntervalUnion Not()const{
                IntervalUnion result;

                /*
                        A) this  \subset world => this
                        B) world \subset this => that

                 */

                if( left.point > right.point )
                        return result;

                if( 0.0 == left.point && left.is_open ){
                        result.children.push_back(Interval::Closed(0.0, 0.0));
                } else if( 0.0 < left.point ){
                        result.children.push_back(Interval{ IntervalEndPoint{false, 0.0}, 
                                                   left.Switch() } );
                }

                if( 1.0 == right.point && right.is_open ){
                        result.children.push_back(Interval::Closed(1.0, 1.0));
                } else if( right.point < 1.0 ){
                        result.children.push_back(Interval{ right.Switch(),       
                                                   IntervalEndPoint{false, 1.0} } );
                }

                return result;
        }

        bool IsSubsetOf(Interval const& that)const{
                if( that.left.point > left.point )
                        return false;
                if( that.left.point ==left.point ){
                        if( that.left.is_open != left.is_open && that.left.is_open )
                                return false;
                }
                if( that.right.point < right.point )
                        return false;
                if( that.right.point ==right.point ){
                        if( that.right.is_open != right.is_open && that.right.is_open )
                                return false;
                }
                return true;
        }

};




struct SortIntervalUnions{
        bool operator()(Union const& a, Union const& b)const{
                if( a.children.size() != b.children.size() )
                        return a.children.size() < b.children.size();
                for(size_t idx=0;idx!=a.children.size();++idx){
                        auto a_interval = boost::get<Interval>(&a.children[idx]);
                        auto b_interval = boost::get<Interval>(&b.children[idx]);
                        BOOST_ASSERT( a_interval );
                        BOOST_ASSERT( b_interval );
                        if( a_interval->left != b_interval->left )
                                return  a_interval->left < b_interval->left;
                        if( a_interval->right != b_interval->right )
                                return  a_interval->right < b_interval->right;
                }
                return false;
        }
};

std::string ToString(BorelSet const& b){
        struct ToStringImpl{
                void operator()(Omega const&){
                        ostr_ << "Omega";
                }
                void operator()(Nul const&){
                        ostr_ << "Nul";
                }
                void operator()(Not const& obj){
                        ostr_ << "Not{";
                        boost::apply_visitor(*this, obj.child);
                        ostr_ << "}";
                }
                void operator()(Union const& obj){
                        ostr_ << "Union{";
                        for(size_t idx=0;idx!=obj.children.size();++idx){
                                if( idx != 0 )
                                        ostr_ << ", ";
                                boost::apply_visitor(*this, obj.children[idx]);
                        }
                        ostr_ << "}";
                }
                void operator()(Intersection const& obj){
                        ostr_ << "Intersection{";
                        for(size_t idx=0;idx!=obj.children.size();++idx){
                                if( idx != 0 )
                                        ostr_ << ", ";
                                boost::apply_visitor(*this, obj.children[idx]);
                        }
                        ostr_ << "}";
                }
                void operator()(Interval const& i){
                        ostr_ << ( i.left.is_open ? "(" : "[" );
                        ostr_ << i.left.point;
                        ostr_ << ",";
                        ostr_ << i.right.point;
                        ostr_ << ( i.right.is_open ? ")" : "]" );
                }
                std::stringstream ostr_;
        };
        ToStringImpl impl;
        boost::apply_visitor(impl, b);
        return impl.ostr_.str();
}

void Display(BorelSet const& b){
        std::cout << ToString(b) << "\n";
}



IntervalUnion ToIntervals(BorelSet const& b){
        struct ToIntervalsImpl : boost::static_visitor<IntervalUnion>{
                IntervalUnion operator()(Omega const&)const{
                        return IntervalUnion{Interval::Closed(0,1)};
                }
                IntervalUnion operator()(Nul const&)const{
                        return IntervalUnion{};
                }
                IntervalUnion operator()(Not const& obj)const{
                        Union result;
                        for(auto const& i : boost::apply_visitor(*this, obj.child).children ){
                                result.Add( i.Not().AsUnion() );
                        }
                        return operator()(result);
                }
                IntervalUnion operator()(Union const& obj)const{
                        IntervalUnion mapped;
                        for(auto const& _ : obj.children ){
                                for( auto const& inner : boost::apply_visitor(*this,_).children ){
                                        mapped.children.push_back(inner);
                                }
                        }
                        std::vector<Interval const*> subs;
                        for(auto const& _ : mapped.children ){
                                subs.push_back(&_);
                        }
                        boost::sort( subs, [](auto const& l, auto const& r){
                                if( l->left.point != r->left.point )
                                        return l->left.point < r->left.point;
                                // prefer closed
                                return l->left.is_open < r->left.is_open;
                        });

                        IntervalUnion result;
                        size_t iter = 0;
                        for(;iter!=subs.size();++iter){
                                if( subs[iter] == 0 )
                                        continue;

                                // for debugging
                                std::vector<Interval const*> dbg_path;
                                dbg_path.push_back(subs[iter]);

                                IntervalEndPoint left  = subs[iter]->left;
                                IntervalEndPoint right = subs[iter]->right;
                                // what is the largest path we can construct
                                subs[iter] = 0;
                                for(size_t j=iter+1;j!=subs.size();){
                                        if( subs[j] == 0 ){
                                                ++j;
                                                continue;
                                        }
                                        auto head = subs[j];

                                        // do these overlap?
                                        if( head->left.point < right.point ||
                                           (head->left.point ==right.point && ! (head->left.is_open && right.is_open) ) ){

                                                // now can we extend right?
                                                if( right.point < head->right.point ||
                                                   (right.point ==head->right.point && right.is_open && !head->right.is_open ) ){
                                                        right = head->right;
                                                        dbg_path.push_back(head);
                                                        // restart nice and slow
                                                        j = iter+1;
                                                        continue;
                                                }
                                        }
                                        ++j;
                                }
                                Interval i{left, right};

                                for(size_t j=iter+1;j!=subs.size();++j){
                                        if( subs[j] == 0 )
                                                continue;
                                        if( subs[j]->IsSubsetOf(i) ){
                                                subs[j] = 0;
                                        }
                                }


                                result.children.push_back(i);

                        }

                        return result;
                }
                IntervalUnion operator()(Intersection const& obj)const{
                        if( obj.children.empty() )
                                return IntervalUnion{};
                        IntervalUnion mapped;
                        for(auto const& _ : obj.children ){
                                for(auto const& inner : boost::apply_visitor(*this,_).children ){
                                        mapped.children.push_back(inner);
                                }
                        }


                        auto first = &mapped.children.at(0);

                        IntervalEndPoint const* upper_left  = &first->left;
                        IntervalEndPoint const* lower_right = &first->right;


                        for(auto const& i : mapped.children ){
                                auto ptr = &i;

                                if( upper_left->point < ptr->left.point )
                                        upper_left = &ptr->left;
                                else if(  upper_left->point == ptr->left.point && 
                                          ! upper_left->is_open && 
                                          ptr->left.is_open )
                                        upper_left = &ptr->left;
                                
                                if( lower_right->point > ptr->right.point )
                                        lower_right = &ptr->right;
                                else if(  lower_right->point == ptr->right.point && 
                                          ! lower_right->is_open && 
                                          ptr->right.is_open )
                                        lower_right = &ptr->right;
                                
                        }

                        std::cout << "*upper_left = " << *upper_left << "\n";
                        std::cout << "*lower_right = " << *lower_right << "\n";

                        if( upper_left->point <  lower_right->point ||
                            ( upper_left->point == lower_right->point 
                              && ! upper_left->is_open && ! lower_right->is_open ) ){
                                return IntervalUnion{ Interval{ *upper_left, *lower_right} };
                        }

                        return IntervalUnion{};
                        
                }
                IntervalUnion operator()(Interval const& i)const{
                        return IntervalUnion{i};
                }
        };

        auto tmp = boost::apply_visitor(ToIntervalsImpl(), b);
        boost::sort(tmp.children);
        return tmp;
}

struct BorelFamily : std::vector<BorelSet>{
        using impl_type = std::vector<BorelSet>;
        template<class... Args>
        BorelFamily(Args&&... args):impl_type{args...}{}
        friend std::ostream& operator<<(std::ostream& ostr, BorelFamily const& self){
                ostr << "{";
                for(size_t idx=0;idx!=self.size();++idx){
                        ostr << ( idx == 0 ? "" : ", " ) << ToString(self[idx]);
                }
                return ostr << "}";
        }
        void Display(std::ostream& out)const{
                for(size_t idx=0;idx!=size();++idx){
                        out << "    " << std::setw(2) << idx << " : " << ToString(at(idx)) << "\n";
                }
        }
};

#if 1
BorelFamily GenerateSigmaAlgebra(BorelFamily const& family){
        BorelFamily head = family;
        std::vector<BorelSet> to_add;


        std::set<IntervalUnion> interval_set;
        for(auto const& _ : family ){
                interval_set.insert( ToIntervals(_) );
        }

        auto test = [&](BorelSet const& b){
                auto iu = ToIntervals(b);
                if( ! interval_set.count( iu ) ){
                        std::cout << "====== found new ======\n";
                        std::cout << "    b  = " << ToString(b) << "\n";
                        std::cout << "    iu = " << ToString(iu.AsUnion()) << "\n";

                        to_add.push_back(b);
                        interval_set.insert( iu );
                        return 1;
                }
                return 0;
        };

        for(;;){
                int changes = 0;
                for(auto const& _ : head ){
                        auto complement = Not{_};
                        changes += test(complement);
                }
                for(size_t i=0;i+1<head.size();++i){
                        for(size_t j=i+1;j<head.size();++j){
                                auto u = Union{ head[i], head[j] };
                                changes += test(u);
                        }
                }
                if( changes == 0 )
                        break;
                for(auto const& _ : to_add )
                        head.push_back(_);
                //break;
        }

        BorelFamily result;
        for(auto const& _ : interval_set ){
                result.push_back( _.AsUnion() );
        }
        return result;



        return head;

}
#endif

int main(){

        BorelSet b = Intersection{ Interval{ Closed(0.0 ), Closed(0.25) },
                            Interval{ Open(0.25), Open(0.50) },
                            Interval{ Closed(0.1), Closed(0.6) } };
        Display(b);


        
        Display(ToIntervals(b).AsUnion());
        
        BorelFamily f0{ Omega(), Nul() };
        std::cout << "f0 = " << f0 << "\n";
        BorelFamily f1{ Omega(), Nul(),
                        Interval{ Closed(0.0), Open(0.5) },
                        Interval{ Closed(0.5), Closed(1.0) } };
        std::cout << "f1 = " << f1 << "\n";
        BorelFamily f2{ Omega(), Nul(),
                        Interval{ Closed(0.0), Open(0.25) },
                        Interval{ Closed(0.25), Open(0.50) },
                        Interval{ Closed(0.50), Open(0.75) },
                        Interval{ Closed(0.75), Closed(1.0) } };
        std::cout << "f2 = " << f2 << "\n";

        std::cout << "GenerateSigmaAlgebra(f2):\n";
        GenerateSigmaAlgebra(f2).Display(std::cout);

        //example_0();
        //example_1();
        //example_2();


}











